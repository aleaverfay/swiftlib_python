# The MIT License (MIT)
#
# Copyright (c) 2014 Andrew Leaver-Fay, Tim Jacobs, Hayretin Yumerefendi, Brian Kuhlman.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
#     in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import blargs
import genetic_code
import amino_acids as aa
import math
import sys

def aastr_for_integer( aaindex ) :
    assert( aaindex >= 0 and aaindex < 21 )
    return aa.amino_acids[aaindex] if aaindex != 20 else "STOP"

# A = 0, C = 1, G = 2, T = 3
# codon = TCG --> 3*16 + 1*4 + 2 = 54
class GeneticCodeMapper :
    def __init__( self ) :
        self.base_to_index = { "A" : 0, "C" : 1, "G" : 2, "T" : 3 }
        self.mapper = [0]*64
        for codon, aastr in genetic_code.genetic_code_codons.iteritems() :
            ci = self.codon_index( codon )
            aaind = 20
            if aastr != "STOP" :
                aastr_up = aastr.upper()
                assert( aastr_up in aa.longer_names )
                aaind = aa.amino_acids.index( aa.longer_names[ aastr_up ] )
            self.mapper[ ci ] = aaind

    def codon_index( self, codon ) :
        index = 0
        for i in xrange( 3 ) :
            index = index*4 + self.base_to_index[ codon[i] ]
        return index

class LexicographicalIterator :
    def __init__( self, dimsizes ) :
        self.size = len(dimsizes)
        self.dimsizes = list( dimsizes )
        self.dimprods = [1] * self.size
        for i in xrange( self.size - 1, 0, -1 ) :
            self.dimprods[ i-1 ] = self.dimprods[ i ] * self.dimsizes[ i ]
        self.search_space_size = self.dimprods[ 0 ] * self.dimsizes[ 0 ]
        self.pos = [0] * self.size
        self.at_end = False

    def increment( self ) :
        i = self.size
        while i > 0 :
            i -= 1
            self.pos[ i ] += 1
            if self.pos[ i ] == self.dimsizes[ i ] :
                self.pos[ i ] = 0
            else :
                return True
        self.at_end = True
        return False

    def upper_diagonal_increment( self ) :
        for i in xrange( self.size - 1, -1, -1 ) :
            self.pos[ i ] += 1
            if ( self.pos[ i ] == self.dimsizes[ i ]  ) :
                self.pos[ i ] = 0
            else :
                beyond_end = False
                for k in xrange(i+1,self.size) :
                    self.pos[ k ] = self.pos[ i ] + k - i
                    if self.pos[k] >= self.dimsizes[k] :
                        beyond_end = True
                        break
                if beyond_end and i == 0 :
                    for k in xrange( self.size ) :
                        self.pos[ k ] = 0
                    self.at_end = True
                    return False
                elif not beyond_end :
                    return True
        self.at_end = True
        return False


    def reset( self ) :
        for i in xrange( self.size ) : self.pos[ i ] = 0
        self.at_end = False

    def upper_diagonal_reset( self ) :
        beyond_end = False
        for i in xrange( self.size ) :
            self.pos[ i ] = i
            if i >= self.dimsizes[i] :
                beyond_end = True
        if beyond_end :
            for i in xrange( self.size ) :
                self.pos[ i ] = 0
            self.at_end = True
        else :
            self.at_end = False

    def index( self ) :
        ''' return the integer index representing the state of the lex'''
        ind = 0
        for i in xrange( self.size ) :
            ind += self.pos[ i ] * self.dimprods[ i ]
        return ind
    def set_from_index( self, ind ) :
        ''' set the state of the lex given a previously computed index'''
        for i in xrange( self.size ):
            self.pos[ i ] = ind / self.dimprods[i]
            ind = ind % self.dimprods[ i ]

class DegenerateCodon :
    def __init__( self ) :
        self.pos = [ [] ] * 3
        for i in xrange(3) : self.pos[i] = [ False ] * 4
        self.which = [ [] ] * 3
        for i in xrange(3) : self.which[i] = []
        self.count_pos = [ 0 ] * 3
        self.degenerate_base_names = {
            (True, False,False,False) : "A",
            (False,True, False,False) : "C",
            (False,False,True, False) : "G",
            (False,False,False,True ) : "T",
            (True, False,False,True ) : "W",
            (False,True, True, False) : "S",
            (True, True, False,False) : "M",
            (False,False,True, True ) : "K",
            (True, False,True, False) : "R",
            (False,True, False,True ) : "Y",
            (False,True, True, True ) : "B",
            (True, False,True, True ) : "D",
            (True, True, False,True ) : "H",
            (True, True, True, False) : "V",
            (True, True, True, True ) : "N" }
        self.names_to_bases = {}
        for bases in self.degenerate_base_names :
            self.names_to_bases[ self.degenerate_base_names[ bases] ] = bases;

    def codon_string( self ) :
        output_codon_string = ""
        for i in xrange(3) :
            output_codon_string += self.degenerate_base_names[ ( self.pos[i][0], self.pos[i][1], self.pos[i][2], self.pos[i][3] ) ]
        return output_codon_string

    def set_pos( self, codon_pos, base ) :
        assert( codon_pos < 3 and codon_pos >= 0 )
        assert( base < 4 and base >= 0 )
        if not self.pos[ codon_pos ][ base ] :
            self.which[ codon_pos ].append( base )
            self.pos[ codon_pos ][ base ] = True
            self.count_pos[ codon_pos ] += 1
    def reset( self ) :
        for i in xrange(3) :
            self.which[i][:] = []
            self.count_pos[i] = 0
            for j in xrange(4) : self.pos[i][j] = False
    def diversity( self ) :
        for i in xrange(3) :
            if self.count_pos[i] == 0 :
                return this.infinity
        div =  1
        for i in xrange(3) :
            div *= self.count_pos[i]
        return div
    def log_diversity( self ) :
        for i in xrange(3) : assert( self.count_pos[i] != 0 )
        ld = 0.0
        for i in xrange(3) : ld += math.log( self.count_pos[i] )
        return ld
    def index_from_lex( self, lex ) :
        """Get the index for a particular codon using a lex that's dimensioned from self.count_pos"""
        codon_index = 0
        for i in xrange(3) :
            codon_index = codon_index * 4 + self.which[i][lex.pos[i]]
        return codon_index
    def set_from_lex( self, lex ) :
        """
        Set the state for this degenerate codon using a lex that's iterating over all (2**4-1)**3 = 3375 codon options.
        Returns False if this is not a reasonable assignment; i.e. not all codon positions contain at least one base.
        """
        self.reset()
        for i in xrange(3) :
            posi = lex.pos[i]+1 # take "14" to mean "all degererate nucleotides" and "0" to mean "only A"
            sigdig = 8
            for j in xrange(4) :
                if posi / sigdig != 0 :
                    self.set_pos( i, 3-j ) # so A = 0 and T = 3
                posi = posi % sigdig
                sigdig /= 2
        return True


class AALibrary :
    def __init__( self ) :
        self.infinity = -1.0;
        self.gcmapper = GeneticCodeMapper()
        self.max_dcs_per_pos = 1
        self.max_oligos_per_stretch = 0
        self.max_oligos_total = 0
        self.n_stretches = 0

    def aas_for_degenerate_codon( self, degenerate_codon ) :
        aas = [ False ] * 21 # 21 because the stop codon counts as a codon.
        lex = LexicographicalIterator( degenerate_codon.count_pos )
        while not lex.at_end :
            codon_index = degenerate_codon.index_from_lex( lex )
            aas[ self.gcmapper.mapper[ codon_index ] ] = True
            lex.increment()
        return aas

    def enumerate_aas_for_all_degenerate_codons( self ) :
        if ( hasattr( self, "aas_for_dc") ) :
            return;
        dims = [ 15, 15, 15 ];
        self.aas_for_dc = [];
        self.diversities_for_dc = []
        self.dclex = LexicographicalIterator( dims )
        dc = DegenerateCodon();
        self.dclex.reset()
        while not self.dclex.at_end :
            dc.set_from_lex( self.dclex );
            self.aas_for_dc.append( self.aas_for_degenerate_codon( dc ) )
            self.diversities_for_dc.append( dc.diversity() )
            self.dclex.increment()


    #format should be a table with N columns and 23 rows
    # row 1 is a header, which just gives the sequence positions
    # row 2 defines stretch boundaries
    # row 3 gives the maximum number of DCs for each position
    # column 1 gives the amino acid names
    # row1/column1 gives nothing
    # all other row/column combinations should be integers
    def load_library( self, libname ) :
        lines = [x.strip() for x in open( libname ).readlines() ]
        assert( len(lines) == 23 )
        row1 = lines[0].split(",")
        self.n_positions = len(row1)-1
        self.aa_counts = [ [] ] * self.n_positions
        self.required  = [ [] ] * self.n_positions
        self.forbidden = [ [] ] * self.n_positions
        self.stretch_reps = self.n_positions * [ 0 ]
        self.max_dcs_for_pos = self.n_positions * [ 1 ]
        for i in xrange(self.n_positions) :
            self.aa_counts[ i ] = [ 0     ] * 21
            self.required[ i ]  = [ False ] * 21
            self.forbidden[ i ] = [ False ] * 21
        self.orig_pos = [ x.strip() for x in row1[1:] ]

        # read out the stretch-boundary data from the CSV input, marking stretch
        # representatives for each position (the first position in a stretch
        # being its own representative) and counting the number of stretches
        row2 = lines[1]
        row2cols = row2.split(",")[1:]
        last_rep = 0

        # the first position is always considered to be the start of a stretch
        self.n_stretches = 1
        self.stretch_reps[0] = 0
        for i in xrange(1,self.n_positions) :
            if row2cols[i] == "|" :
                last_rep = i
                self.n_stretches += 1
            self.stretch_reps[i] = last_rep

        # increase the maximum number of oligos total if it is less than the
        # total number of stretches, othrwise, there wouldn't be enough oligos
        # to cover all the stretches and the optimization would be nonsensical
        if self.max_oligos_total < self.n_stretches :
            self.max_oligos_total = self.n_stretches

        # read out the per-position counts of the number of degenerate codons
        # allowed at each position from the CSV input text
        row3 = lines[2]
        row3cols = row3.split(",")[1:]
        self.max_dcs_per_pos = 1
        for i in xrange( self.n_positions ) :
            self.max_dcs_for_pos[ i ] = int( row3cols[i] )
            if self.max_dcs_for_pos[i] > self.max_dcs_per_pos :
                self.max_dcs_per_pos = self.max_dcs_for_pos[i]

        # increase the maximum number of oligos per stretch to give space for
        # at least one position in a stretch to have as many degenerate codons
        # as has been allowed at that position, otherwise, there is no point in
        # entertaining that many degenerate codons at all
        if self.max_oligos_per_stretch == 0 :
            self.max_oligos_per_stretch = self.max_oligos_total - self.n_stretches + 1
        elif self.max_oligos_per_stretch < self.max_dcs_per_pos :
            self.max_oligos_per_stretch = self.max_dcs_per_pos

        # increase the maximum number of oligos total to give space for the maximum number
        # of oligos in a single stretch - 1 + the total number of stretches
        # otherwise, there is no point in looking at multiple degenerate codons at a single
        # position or multiple oligos to cover a single stretch
        if self.max_oligos_total < self.n_stretches + self.max_oligos_per_stretch - 1 :
            self.max_oligos_total = self.n_stretches + self.max_oligos_per_stretch - 1

        self.max_per_position_error = 0
        obs = [ 0 ] * self.n_positions
        for i in xrange(20) :
            line = lines[ i + 3 ]
            vals = line.split(",")[1:]
            for j in xrange(len(vals)):
                if vals[j] == "!" :
                    self.forbidden[j][i] = True
                elif vals[j] == "*" :
                    self.required[j][i] = True
                else :
                    self.aa_counts[j][i] = int(vals[j])
                obs[ j ] += self.aa_counts[j][i]
        for i in xrange( self.n_positions ) :
            if obs[i] > self.max_per_position_error :
                self.max_per_position_error = obs[i]
    def error_given_aas_for_pos( self, pos, aas ) :
        error = 0
        for i in xrange(20) :
            icount = self.aa_counts[ pos ][ i ]
            if not aas[ i ] :
                if self.required[ pos ][ i ] :
                    return self.infinity
                else :
                    error += icount
            else :
                if self.forbidden[ pos ][ i ] :
                    return self.infinity
                if icount < 0 :
                    error -= icount
        return error

    def error_given_aas_for_pos_ignore_req( self, pos, aas ) :
        error = 0
        for i in xrange(20) :
            icount = self.aa_counts[ pos ][ i ]
            if not aas[ i ] :
                if icount > 0 :
                    error += icount
            else :
                if self.forbidden[ pos ][ i ] :
                    return self.infinity
                if icount < 0 :
                    error -= icount
        return error

    def useful_aaind_for_pos( self, aas, pos ) :
        aaind = 0
        for i in xrange(21) :
            iuseful = self.aa_counts[pos][i] != 0 or self.required[pos][i];
            aaind = 2*aaind + ( 1 if aas[i] and iuseful else 0 )
        return aaind

    def find_useful_codons( self ) :
        self.useful_codons = []

        # the diversities for the degenerate codons; this is a three-dimensional array
        # index 0: which position
        # index 1: what error
        # index 2: either 0 (for the codon's diveristy) or 1 (for the codon's index)

        div_for_codons = []
        for i in xrange( self.n_positions ) :
            self.useful_codons.append( [] )
            div_for_codons.append( {} )
        for i in xrange( 3375 ) :
            iaas = self.aas_for_dc[ i ]
            idiv = self.diversities_for_dc[ i ]
            for j in xrange( self.n_positions ) :
                ijerror = self.error_given_aas_for_pos_ignore_req( j, iaas )
                if ( ijerror == self.infinity ) : continue;
                ij_aaind = self.useful_aaind_for_pos( iaas, j )

                # keep this codon if it's the first codon producing ijerror
                # OR it's the first codon with the given aa index producing ijerror
                # OR the diversity for this codon is smaller than the smallest-seen
                # diversity of any codon producing ijerror with the given aa index.
                if ( ijerror not in div_for_codons[j] or
                     ij_aaind not in div_for_codons[j][ ijerror ] or
                     div_for_codons[j][ ijerror][ ij_aaind ][ 0 ] > idiv ) :
                    if ( ijerror not in div_for_codons[j] )  :
                        div_for_codons[j][ ijerror ] = {}
                    div_for_codons[j][ ijerror ][ ij_aaind ] = [ idiv, i ]; # JS->Py note: could be a tuple?
        for i in xrange( self.n_positions ) :
            for j in div_for_codons[ i ] :
                if ( j not in div_for_codons[i] ) : continue # JS->PY note: maybe unnecessary?
                jaainds = div_for_codons[i][j]
                for k in jaainds :
                    self.useful_codons[i].append( jaainds[k][1] )


    def codon_inds_from_useful_codon_lex( self, position, useful_codon_lex ) :
        inds = []
        for i in xrange( len( useful_codon_lex.pos ) ) :
            inds.append( self.useful_codons[ position ][ useful_codon_lex.pos[i] ] )
        return inds

    def compute_smallest_diversity_for_all_errors( self ) :

        self.enumerate_aas_for_all_degenerate_codons();
        self.find_useful_codons();

        self.divmin_for_error_for_n_dcs = [ [] ] * self.n_positions
        self.codons_for_error_for_n_dcs = [ [] ] * self.n_positions
        self.errors_for_n_dcs_for_position = [ [] ] * self.n_positions
        for i in xrange(self.n_positions) :
            #print "compute_smallest_diversity_for_all_errors", i, len( self.divmin_for_error_for_n_dcs ), len( self.max_dcs_for_pos )
            self.divmin_for_error_for_n_dcs[i] =    ( self.max_dcs_for_pos[i]+1 ) * [ [] ]
            self.codons_for_error_for_n_dcs[i] =    ( self.max_dcs_for_pos[i]+1 ) * [ [] ]
            self.errors_for_n_dcs_for_position[i] = ( self.max_dcs_for_pos[i]+1 ) * [ [] ]
            for j in xrange(1, self.max_dcs_for_pos[i]+1) :
                self.divmin_for_error_for_n_dcs[i][j] = ( self.max_per_position_error+1 ) * [ self.infinity ]
                self.codons_for_error_for_n_dcs[i][j] = ( self.max_per_position_error+1 ) * [ self.infinity ]
                self.errors_for_n_dcs_for_position[i][j] = []
        aas_for_combo = 21 * [ False ]
        for i in xrange( self.n_positions ) :
            for j in xrange(1,self.max_dcs_for_pos[i]+1) :
                dims = j*[ 0 ]
                print "finding smallest diversity codons for errors", i, j, ":"
                for k in xrange(j) :
                    dims[k] = len( self.useful_codons[i] )
                    print dims[k],
                print
                jlex = LexicographicalIterator( dims )
                jlex.upper_diagonal_reset()

                while not jlex.at_end :
                    for k in xrange(21) : aas_for_combo[k] = False
                    diversity = 0
                    for k in xrange(j) :
                        kcodon = self.useful_codons[i][ jlex.pos[k] ]
                        diversity += self.diversities_for_dc[kcodon]
                        kaas = self.aas_for_dc[ kcodon ]
                        for l in xrange(21) :
                            aas_for_combo[l] = aas_for_combo[l] or kaas[l]
                    log_diversity = math.log( diversity )
                    error = self.error_given_aas_for_pos( i, aas_for_combo )
                    if ( error != self.infinity ) :
                        #print i, j, error, self.max_per_position_error
                        sys.stdout.flush()
                        prev_diversity = self.divmin_for_error_for_n_dcs[i][j][error]
                        if ( prev_diversity == self.infinity ) :
                            self.errors_for_n_dcs_for_position[i][j].append( error )
                        if ( prev_diversity == self.infinity or prev_diversity > log_diversity ) :
                            self.divmin_for_error_for_n_dcs[i][j][error] = log_diversity
                            self.codons_for_error_for_n_dcs[i][j][error] = self.codon_inds_from_useful_codon_lex( i, jlex )
                    jlex.upper_diagonal_increment();
                self.errors_for_n_dcs_for_position[i][j].sort()

        # self.divmin_for_error = [ [] ] * self.n_positions
        # for i in xrange( self.n_positions ) : self.divmin_for_error[i] = [ ( self.infinity, 0 ) ] * self.max_per_position_error
        # dims = [ 2**4 ] * 3 # i.e. [ 16 ] * 3
        # self.dclex = LexicographicalIterator( dims )
        # dc = DegenerateCodon()
        # for i in xrange( self.n_positions ) :
        #     self.dclex.reset()
        #     while not self.dclex.at_end :
        #         if dc.set_from_lex( self.dclex ) :
        #             aas = self.aas_for_degenerate_codon( dc )
        #             error = self.error_given_aas_for_pos( i, aas )
        #             log_diversity = dc.log_diversity()
        #             prev_diversity = self.divmin_for_error[ i ][ error ][0]
        #             if prev_diversity == self.infinity or log_diversity < prev_diversity :
        #                 # store the diversity and information on the degenerate codon that
        #                 # produced this level of error
        #                 self.divmin_for_error[i][ error ] = ( log_diversity, self.dclex.index() )
        #         self.dclex.increment()


    def optimize_library( self, diversity_cap ) :
        '''
        Run a dynamic programming algorithm to determine the minimum diversity for
        every error level, and return an array of error levels for each position
        that describes the library that fits under the diversity cap with the smallest
        error.  This array can be used with the previously-computed divmin_for_error
        array to figure out which codons should be used at every position.
        '''

        assert( hasattr( self, 'divmin_for_error' ) )
        self.dp_divmin_for_error = [[]] * self.n_positions

        # dp_traceback is an array of tuples:
        #  pos0 = which error level ( from self.divmin_for_error ) for this position
        #  pos1 = which error total ( from self.dp_divmin_for_error ) for the previous position
        self.dp_traceback = [[]] * self.n_positions
        self.error_span = self.max_per_position_error * self.n_positions
        for i in xrange( self.n_positions ) :
            self.dp_divmin_for_error[i] = [ self.infinity ] * self.error_span
            self.dp_traceback[i] = [ ( self.infinity, self.infinity ) ] * self.error_span

        # take care of position 0: copy self.divmin_for_error[0] into self.dp_divmin_for_eror
        for i in xrange( self.max_per_position_error ) :
            self.dp_divmin_for_error[0][i] = self.divmin_for_error[0][i][0]
            self.dp_traceback[0][i] = ( i, 0 ) # traceback doesn't proceed beyond position 0

        for i in xrange( 1, self.n_positions )  :
            # solve the dynamic programming problem for residues 0..i
            for j in xrange( self.error_span ) :
                j_divmin = self.infinity
                j_traceback = None
                for k in xrange( min( j, self.max_per_position_error ) ) :
                    if self.dp_divmin_for_error[i-1][j-k] == self.infinity : continue
                    if self.divmin_for_error[i][k][0] == self.infinity : continue
                    divsum = self.divmin_for_error[i][k][0] + self.dp_divmin_for_error[i-1][j-k]
                    if j_divmin == self.infinity or divsum < j_divmin :
                        j_divmin = divsum
                        j_traceback = ( k, j-k )
                if j_divmin != self.infinity :
                    self.dp_divmin_for_error[i][j] = j_divmin
                    self.dp_traceback[i][j] = j_traceback
        return self.traceback()

    def traceback( self ) :
        # now the traceback
        optimal_error_traceback = [ 0 ] * self.n_positions
        log_diversity_cap =  math.log( diversity_cap )

        #for i in xrange( self.error_span ) :
        #    if self.dp_divmin_for_error[-1][i] != self.infinity :
        #        print "Error of",i,"requires diversity of %5.3f" % self.dp_divmin_for_error[-1][i]

        for i in xrange( self.error_span ) :
            if self.dp_divmin_for_error[-1][i] != self.infinity and self.dp_divmin_for_error[-1][i] < log_diversity_cap :
                best = i
                print "Minimum error of", i, "with log(diversity) of",self.dp_divmin_for_error[-1][i]
                break
        return self.traceback_from_error_level( best, optimal_error_traceback )

    def traceback_from_error_level( self, error_level, error_traceback ) :
        position_error = [0] * self.n_positions
        for i in xrange( self.n_positions - 1, -1, -1 ) :
            tb = self.dp_traceback[ i ][ error_level ];
            error_traceback[ i ] = tb[0]
            error_level = tb[1]
            position_error[i] = tb[0]
        for i in xrange( self.n_positions ) :
            print "Traceback position", i, "minimum error=", position_error[i]

        return error_traceback

def final_codon_string( position, degenerate_codon, library ) :
    # three things we need:
    # 1: the codon
    # 2: the amino acids that are represented
    # 2b: the counts from the original set of observations for each of the represented aas
    # 3: the amino acids and their counts in the original set of observations that are not represented
    aas_present  = library.aas_for_degenerate_codon( degenerate_codon )
    orig_obs = library.aa_counts[ position ]

    orig_pos_string = "Position %4s" % library.orig_pos[ position ]

    codon_string = ""
    for i in xrange(3) :
        idcpos = degenerate_codon.pos[i]
        base_tuple = ( idcpos[0], idcpos[1], idcpos[2], idcpos[3] )
        codon_string += degenerate_codon.degenerate_base_names[ base_tuple ]

    present_string = ""
    for i in xrange(len(aas_present)) :
        if aas_present[i] :
            present_string += " " + aastr_for_integer( i )
            if i != 20 : present_string += "(" + str(orig_obs[i]) + ")"
    absent_string = "Absent"
    for i in xrange(len(orig_obs)) :
        if orig_obs[i] != 0 and not aas_present[i]:
            absent_string += " " + aastr_for_integer( i ) + "(" + str(orig_obs[i]) + ")"
    log_diversity_string = "log(diversity)= %5.3f" % ( degenerate_codon.log_diversity() )

    return orig_pos_string + " : " + codon_string + " : " + log_diversity_string + " : " + present_string + " : " + absent_string

def print_output_codons( library, error_sequence ) :
    dc = DegenerateCodon()
    diversity_sum = 0
    for i in xrange(library.n_positions) :
        lexind = library.divmin_for_error[ i ][ error_sequence[ i ] ][ 1 ]
        library.dclex.set_from_index(lexind)
        dc.set_from_lex( library.dclex )
        diversity_sum += dc.log_diversity()
        print final_codon_string( i, dc, library )
    print "Max log diversity: ", math.log( diversity_cap ), "Theorical diversity", diversity_sum

def practice_code( library ) :
    print library.max_per_position_error
    print library.aa_counts

    dims = [ 2, 1, 1 ]
    lex = LexicographicalIterator( dims )
    print lex.pos
    lex.increment()
    print lex.pos

    print "bigger lex"
    dims = [ 4, 3, 5 ]
    lex = LexicographicalIterator( dims )
    print lex.pos, lex.dimprods, lex.search_space_size
    for i in xrange(50) : lex.increment()
    print lex.pos
    index = lex.index()
    lex.reset()
    print index, lex.pos
    lex.set_from_index( index )
    print lex.pos

    ttg_degenerate = DegenerateCodon()
    ttg_degenerate.set_pos(0,3); ttg_degenerate.set_pos(1,3); ttg_degenerate.set_pos(2,2);
    aas_ttg = library.aas_for_degenerate_codon( ttg_degenerate )
    for i in xrange( 21 ) :
        if aas_ttg[i] :
            print aastr_for_integer( i )

    ttg_degenerate.set_pos(0,1)
    aas_ttg = library.aas_for_degenerate_codon( ttg_degenerate )
    for i in xrange( 21 ) :
        if aas_ttg[i] :
            print aastr_for_integer( i )

    #ttg_degenerate[1][0] = True; # now try C/T A/T G
    ttg_degenerate.set_pos(1,0)
    aas_ttg = library.aas_for_degenerate_codon( ttg_degenerate )
    for i in xrange( 21 ) :
        if aas_ttg[i] :
            print aastr_for_integer( i )


if __name__ == "__main__" :
    with blargs.Parser( locals() ) as p :
        p.str( "input_csv" ).required()
        p.float( "diversity_cap" ).required()

    print "Loading library"
    library = AALibrary()
    library.load_library( input_csv )

    print "Computing minimum diversity for each error level for each position"
    library.compute_smallest_diversity_for_all_errors()

    print "Running dynamic programming to minimize error while coming under the diversity cap"
    optimal = library.optimize_library( diversity_cap )

    print_output_codons( library, optimal )
