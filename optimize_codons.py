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

    def reset( self ) :
        for i in xrange( self.size ) : self.pos[ i ] = 0
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
            (True, False,True, False) : "P",
            (False,True, False,True ) : "Y",
            (False,True, True, True ) : "B",
            (True, False,True, True ) : "D",
            (True, True, False,True ) : "H",
            (True, True, True, False) : "V",
            (True, True, True, True ) : "N" }
        

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
        Set the state for this degenerate codon using a lex that's iterating over all 2**4**3 = 4096 codon options.
        Returns False if this is not a reasonable assignment; i.e. not all codon positions contain at least one base.
        """
        for i in xrange(3) :
            if lex.pos[i] == 0 : return False

        self.reset()
        for i in xrange(3) :
            posi = lex.pos[i]
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

    def aas_for_degenerate_codon( self, degenerate_codon ) :
        aas = [ False ] * 21 # 21 because the stop codon counts as a codon.
        lex = LexicographicalIterator( degenerate_codon.count_pos )
        while not lex.at_end :
            codon_index = degenerate_codon.index_from_lex( lex )
            aas[ self.gcmapper.mapper[ codon_index ] ] = True
            lex.increment()
        return aas

    #format should be a table with N columns and 21 rows
    # row 1 is a header, which just gives the sequence positions
    # column 1 gives the amino acid names
    # row1/column1 gives nothing
    # all other row/column combinations should be integers
    def load_library( self, libname ) :
        lines = open( libname ).readlines()
        assert( len(lines) == 21 )
        row1 = lines[0].split(",")
        self.n_positions = len(row1)-1
        self.aa_counts = [ [] ]*(self.n_positions)
        for i in xrange(self.n_positions) : self.aa_counts[ i ] = [0] * 20
        self.orig_pos = [ x.strip() for x in row1[1:] ]
        self.max_obs = 0
        for i in xrange(20) :
            line = lines[ i + 1 ]
            vals = line.split(",")[1:]
            iiobs = 0
            for j in xrange(len(vals)):
                self.aa_counts[j][i] = int(vals[j])
                iiobs += self.aa_counts[j][i]
            if iiobs > self.max_obs :
                self.max_obs = iiobs
    def error_given_aas_for_pos( self, pos, aas ) :
        error = 0
        for i in xrange(20) :
            if not aas[ i ] :
                error += self.aa_counts[ pos ][ i ]
        return error

    def compute_smallest_diversity_for_all_errors( self ) :
        self.divmin_for_error = [ [] ] * self.n_positions
        for i in xrange( self.n_positions ) : self.divmin_for_error[i] = [ ( self.infinity, 0 ) ] * self.max_obs
        dims = [ 2**4 ] * 3 # i.e. [ 16 ] * 3
        self.gclex = LexicographicalIterator( dims )
        gc = DegenerateCodon()
        for i in xrange( self.n_positions ) :
            self.gclex.reset()
            while not self.gclex.at_end :
                if gc.set_from_lex( self.gclex ) :
                    aas = self.aas_for_degenerate_codon( gc )
                    error = self.error_given_aas_for_pos( i, aas )
                    log_diversity = gc.log_diversity()
                    prev_diversity = self.divmin_for_error[ i ][ error ][0]
                    if prev_diversity == self.infinity or log_diversity < prev_diversity :
                        # store the diversity and information on the degenerate codon that
                        # produced this level of error
                        self.divmin_for_error[i][ error ] = ( log_diversity, self.gclex.index() )
                self.gclex.increment()

                
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
        self.error_span = self.max_obs * self.n_positions
        for i in xrange( self.n_positions ) :
            self.dp_divmin_for_error[i] = [ self.infinity ] * self.error_span
            self.dp_traceback[i] = [ ( self.infinity, self.infinity ) ] * self.error_span

        # take care of position 0: copy self.divmin_for_error[0] into self.dp_divmin_for_eror
        for i in xrange( self.max_obs ) :
            self.dp_divmin_for_error[0][i] = self.divmin_for_error[0][i][0]
            self.dp_traceback[0][i] = ( i, 0 ) # traceback doesn't proceed beyond position 0

        for i in xrange( 1, self.n_positions )  :
            # solve the dynamic programming problem for residues 0..i
            for j in xrange( self.error_span ) :
                j_divmin = self.infinity
                j_traceback = None
                for k in xrange( min( j, self.max_obs ) ) :
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
        igcpos = degenerate_codon.pos[i]
        base_tuple = ( igcpos[0], igcpos[1], igcpos[2], igcpos[3] )
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
    gc = DegenerateCodon()
    diversity_sum = 0
    for i in xrange(library.n_positions) :
        lexind = library.divmin_for_error[ i ][ error_sequence[ i ] ][ 1 ]
        library.gclex.set_from_index(lexind)
        gc.set_from_lex( library.gclex )
        diversity_sum += gc.log_diversity()
        print final_codon_string( i, gc, library )
    print "Max log diversity: ", math.log( diversity_cap ), "Theorical diversity", diversity_sum

def practice_code( library ) :
    print library.max_obs
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
