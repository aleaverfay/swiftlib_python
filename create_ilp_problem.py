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
import math
import optimize_codons

class DegenerateCodonForILP :
   def __init__( self ) :
      self.n_dcs = 1
      self.dcs = []
      self.stretch = 1
      self.position = 1
      self.index = 1
      self.index_by_ndcs = 1
      self.name = ""
      self.log_div = 0
      self.error = 0
      self.codons = ""
      self.aas = ""

   def __str__( self ) :
      s = "DC: " + self.name + " on " + str(self.position) + " of stretch: " + str(self.stretch) \
          + " ndcs: " + str( self.n_dcs ) + " ind: " + str( self.index ) + " div: " + str( self.log_div ) \
          + " err: " + str( self.error )
      return s

class NDCCountVar :
   def __init__( self ) :
      self.name = ""
      self.stretch = 0
      self.n_dcs = 0
      self.pos = 0
      self.aux_vars_part_of = []
   def __str__( self ) :
      s = "NDCCountVar " + self.name + " " + str(self.stretch) + " part of:\n"
      for avpo in self.aux_vars_part_of :
         s += "  " + avpo + "\n"
      return s
      
class StretchNPrimersVar :
   def __init__( self ) :
      self.name = ""
      self.index = 0
      self.nprimers = 0
      self.aux_vars = []

class InequalityAuxVar :
   def __init__( self ) :
      self.name = ""
      self.c0 = 0

def define_variables_for_library( library ) :
   vars = {}

   print "Define variables for library, max_dcs_per_pos", library.max_dcs_per_pos
   if ( library.max_dcs_per_pos == 1 ) : return vars

   ndcs_per_pos = []
   count_stretch = -1
   for i in xrange( library.n_positions ) :
      if library.stretch_reps[i] == i :
         count_stretch += 1
      posi_ndc_vars = []
      for j in xrange( 1, library.max_dcs_for_pos[i]+1 ) :
         jvar = NDCCountVar()
         jvar.name = "NDCS_POS%d_IS_%d" % ( i, j )
         jvar.stretch = count_stretch
         jvar.n_dcs = j
         jvar.pos = i
         posi_ndc_vars.append( jvar )
      ndcs_per_pos.append( posi_ndc_vars )
   vars[ "ndcs" ] = ndcs_per_pos

   stretch_members = []
   curr_stretch_pos = []
   count_stretch = -1
   for i in xrange( library.n_positions ) :
      if library.stretch_reps[i] == i :
         count_stretch += 1
         if count_stretch != 0 :
            stretch_members.append( curr_stretch_pos )
            curr_stretch_pos = []
      curr_stretch_pos.append( i )
   stretch_members.append( curr_stretch_pos )

   stretch_combos = []
   all_gt_aux_vars = []
   all_lt_aux_vars = []
   for i in xrange( count_stretch+1 ) :
      curr_combos = []
      imembs = len( stretch_members[i] )
      dims = imembs * [ 0 ]
      for j in xrange( imembs ) :
         jpos = stretch_members[i][j]
         dims[j] = library.max_dcs_for_pos[ jpos ]
      lex = optimize_codons.LexicographicalIterator( dims )
      while not lex.at_end :
         nprimers_var = StretchNPrimersVar()
         varname = "S%d_COMBO" % i
         size = 1
         for j in xrange( imembs ) :
            jsize = lex.pos[j] + 1
            varname += "_%d" % ( jsize )
            size *= jsize
         nprimers_var.name = varname
         nprimers_var.nprimers = size

         gt_auxvar = InequalityAuxVar()
         gt_auxvar.name = "AND_GT_" + varname
         gt_auxvar.c0 = imembs - 1
         gt_auxvars = [ gt_auxvar ]
         lt_auxvars = []
         nprimers_var.aux_vars.append( gt_auxvars[0].name )
         for j in xrange( imembs ) :
            jpos = stretch_members[i][j]
            lt_auxvar = InequalityAuxVar()
            lt_auxvar.name = ( "AND_LT_%d_" % jpos ) + varname
            lt_auxvar.c0 = 0
            lt_auxvars.append( lt_auxvar )
            ndcs_per_pos[ jpos ][ lex.pos[j] ].aux_vars_part_of.append( gt_auxvar.name )
            ndcs_per_pos[ jpos ][ lex.pos[j] ].aux_vars_part_of.append( lt_auxvars[-1].name )
            nprimers_var.aux_vars.append( lt_auxvars[-1].name )
        
         all_gt_aux_vars.append( gt_auxvar )
         all_lt_aux_vars.extend( lt_auxvars )

         curr_combos.append( nprimers_var )
         lex.increment()
      stretch_combos.append( curr_combos )
   vars[ "stretch_combos" ] = stretch_combos
   vars[ "gt_aux_vars" ] = all_gt_aux_vars
   vars[ "lt_aux_vars" ] = all_lt_aux_vars

   return vars
      

def create_dcs_for_ilp_for_pos( library, pos, stretch ) :
   degcodon = optimize_codons.DegenerateCodon()
   dclist = []
   ncodons_by_ndc = []
   count = 0
   for i in xrange( 1, library.max_dcs_for_pos[ pos ]+1 ) :
      count_by_ndcs = 0
      posi_divmin = library.divmin_for_error_for_n_dcs[pos][i]
      print i, library.max_dcs_for_pos[ pos ], len( posi_divmin ), library.max_per_position_error
      for j in xrange( library.max_per_position_error + 1 ) :
         if posi_divmin[ j ] != library.infinity :
            count += 1
            count_by_ndcs += 1
            dc = DegenerateCodonForILP()
            dc.name = "DC_" + str(pos) + "_"  + str(count)
            dc.stretch = stretch
            dc.position = pos
            dc.index = count
            dc.index_by_ndcs = count_by_ndcs
            dc.n_dcs = i
            dc.log_div = posi_divmin[ j ]
            dc.error = j
            dclist.append( dc )
            codon_inds = library.codons_for_error_for_n_dcs[ pos ][ i ][ j ]
            jcodons = []
            aas_present = 21 * [ False ]
            for k in codon_inds :
               library.dclex.set_from_index( k )
               degcodon.set_from_lex( library.dclex )
               jcodons.append( degcodon.codon_string() )
               kaas = library.aas_for_dc[ k ]
               for l in xrange(21) :
                  aas_present[ l ] |= kaas[ l ]
            dc.codons = "+".join( jcodons )
            present_aas = []
            for k in xrange(21) :
               if aas_present[k] :
                  present_aas.append( optimize_codons.aastr_for_integer( k ))
            dc.aas = "".join( present_aas )
      #print "Appending to ncodons_by_ncd: ", count_by_ndcs, " ncodons_by_ndc: ", ncodons_by_ndc
      ncodons_by_ndc.append( count_by_ndcs )
   # for dc in dclist :
   #    print dc
   # print
   # print "Ncodons by ndc: ", ncodons_by_ndc
   return dclist, ncodons_by_ndc


def dc_column( dc, mdcs ) :
   line1 = " " + dc.name + " error " + str( dc.error ) + "\n"
   line2 = " " + dc.name + " libsize " + str( "%6f" % dc.log_div ) + "\n"
   line3 = " " + dc.name + ( " onlyone_%d" % dc.position ) + " 1\n"
   if ( mdcs ) :
      line4 = " " + dc.name + ( " ndc_%d_%d_or_everything" % ( dc.position, dc.n_dcs ) ) + " 1\n"
      line5 = " " + dc.name + ( " ndc_%d_%d_or_%d" % ( dc.position, dc.n_dcs, dc.index_by_ndcs ) ) + " 1\n"
      return [ line1, line2, line3, line4, line5 ]
   else :
      return [ line1, line2, line3 ]

def ndc_count_var_column( ndc_var, ndcs  ) :
   #print "ndc_count_var_column", ndc_var.name, " ", ndcs, ( "%d\n" % (-1*ndcs) )
   lines = [ " " + ndc_var.name + ( " ndc_%d_%d_or_everything -1\n" % ( ndc_var.pos, ndc_var.n_dcs ) ) ]
   for i in xrange(1,ndcs+1) :
      lines.append( " " + ndc_var.name + ( " ndc_%d_%d_or_%d" % ( ndc_var.pos, ndc_var.n_dcs, i ) ) + " -1\n" )
   # these aux vars represent the boolean AND logic when trying to compute
   # the product of the number of DCs used at all positions in a single stretch
   # which is needed to get the number of primers demanded by the stretch
   for aux in ndc_var.aux_vars_part_of :
      lines.append( " " + ndc_var.name + " " + aux + " 1\n" )
   return lines

def stretch_n_primers_column( snpv ) :
   lines = [ " " + snpv.name + ( " total_num_primers %d\n" % snpv.nprimers ) ]
   for aux in snpv.aux_vars :
      lines.append( " " + snpv.name + " " + aux + " -1\n" )
   return lines

# def columns_for_position( library, pos, dc ) :
#    lines = []
#    varlist = []
#    count = 0
#    divmin_for_error_for_pos = library.divmin_for_error[ pos ]
#    for i in xrange( library.max_obs ) :
#        if divmin_for_error_for_pos[ i ][ 0 ] != library.infinity :
#            count += 1
#            varname = "DC_" + str(count) + "_AT_" + str(pos)
#            line1 = " " + varname + " error " + str(i) + "\n"
#            line2 = " " + varname + " libsize " + ( "%6f" % divmin_for_error_for_pos[ i ][ 0 ] ) + "\n"
#            line3 = " " + varname + " onlyone_" + str(pos) + " 1\n"
#            lines.append( line1 )
#            lines.append( line2 )
#            lines.append( line3 )
#            library.dclex.set_from_index( divmin_for_error_for_pos[ i ][ 1 ] )
#            dc.set_from_lex( library.dclex )
#            varlist.append( ( varname, dc.codon_string() ))
#            #varmap.append( varname + " " +  dc.codon_string() + "\n" )
#    return varlist, lines


# python create_ilp_problem.py --input_csv the_problem.csv --diversity_cap 1e7

if __name__ == "__main__" :
   with blargs.Parser( locals() ) as p :
       p.str( "input_csv" ).required()
       p.float( "diversity_cap" ).required()
       p.int( "nprimer_limit" ).default(-1)
       p.str( "ilp_input_prefix" ).required()

   print "Loading library"
   library = optimize_codons.AALibrary()
   library.load_library( input_csv )

   print "Finding minimum diversity codons"
   library.compute_smallest_diversity_for_all_errors()

   count_stretch = -1
   dcs = []
   ncodons_by_ndcs = []
   for i in xrange( library.n_positions ) :
      if library.stretch_reps[i] == i :
         count_stretch += 1
      idcs, ncodons_by_ndc = create_dcs_for_ilp_for_pos( library, i, count_stretch )
      dcs.append( idcs )
      ncodons_by_ndcs.append( ncodons_by_ndc )

   if ( nprimer_limit == -1 ) : nprimer_limit = count_stretch;

   #print "NCODONS_BY_NDCS: ", ncodons_by_ndcs

   vars = define_variables_for_library( library )

   print "beginning to write ILP problem"
   fmps_lines = []
   fmps_lines.append( "NAME MULTIPLE_DEGENERATE_CODON_OPTIMIZATION\n" )
   fmps_lines.append( "ROWS\n" )
   fmps_lines.append( " N error\n" )
   fmps_lines.append( " L libsize\n" )

   for i in xrange( library.n_positions ) :
      fmps_lines.append( " E onlyone_" + str(i) + "\n" )

   if ( library.max_dcs_per_pos != 1 ) :
       fmps_lines.append( " L total_num_primers\n" )
       for aux in vars[ "gt_aux_vars" ] :
          fmps_lines.append( " L " + aux.name + "\n" )
       for aux in vars[ "lt_aux_vars" ] :
          fmps_lines.append( " G " + aux.name + "\n" )
       for i in xrange( library.n_positions ) :
          for j in xrange( 1, library.max_dcs_for_pos[i] + 1 ) :
             fmps_lines.append( " G ndc_" + str(i) + "_" + str(j) + "_or_everything\n" )
       for i in xrange( library.n_positions ) :
          for j in xrange( 1, library.max_dcs_for_pos[i] + 1 ) :
             for k in xrange( 1, ncodons_by_ndcs[ i ][ j-1 ] + 1 ) :
                fmps_lines.append( " L ndc_%d_%d_or_%d\n" % ( i, j, k ))


   fmps_lines.append( "COLUMNS\n" )
   for i in xrange( library.n_positions ) :
      for dc in dcs[ i ] :
         fmps_lines.extend( dc_column( dc, library.max_dcs_per_pos != 1 ) )

   if ( library.max_dcs_per_pos != 1 ) :
       ndc_vars = vars[ "ndcs" ]
       for i in xrange( library.n_positions ) :
          for j in xrange( library.max_dcs_for_pos[ i ] ) :
             #print "WHAT?!", ncodons_by_ndcs, ncodons_by_ndcs[i], ncodons_by_ndcs[i][j]
             fmps_lines.extend( ndc_count_var_column( ndc_vars[i][j], ncodons_by_ndcs[i][j] ))
    
       stretch_combos = vars[ "stretch_combos" ]
       for sclist in stretch_combos :
          for sc in sclist :
             fmps_lines.extend( stretch_n_primers_column( sc ) )

   fmps_lines.append( "RHS\n" )
   fmps_lines.append( " RHS1 libsize " + ( "%6f\n" % math.log( diversity_cap )))

   if ( library.max_dcs_per_pos != 1 ) :
      fmps_lines.append( " RHS1 total_num_primers " + str(nprimer_limit) + "\n" )
   for i in xrange( library.n_positions ) :
      fmps_lines.append( " RHS1 onlyone_%d 1\n" % i )

   if ( library.max_dcs_per_pos != 1 ) :
       for aux in vars[ "gt_aux_vars" ] :
          fmps_lines.append( " RHS1 " + aux.name + " " + str( aux.c0 ) + "\n" )
       for aux in vars[ "lt_aux_vars" ] :
          fmps_lines.append( " RHS1 " + aux.name + " " + str( aux.c0 ) + "\n" )
       for i in xrange( library.n_positions ) :
          for j in xrange( 1, library.max_dcs_for_pos[i] + 1 ) :
             fmps_lines.append( " RHS1 ndc_" + str(i) + "_" + str(j) + "_or_everything 0\n" )
       for i in xrange( library.n_positions ) :
          for j in xrange( 1, library.max_dcs_for_pos[i] + 1 ) :
             for k in xrange( 1, ncodons_by_ndcs[ i ][ j-1 ] + 1 ) :
                fmps_lines.append( " RHS1 ndc_%d_%d_or_%d 0 \n" % ( i, j, k ))

   fmps_lines.append( "BOUNDS\n" )
   for i in xrange( library.n_positions ) :
      for dc in dcs[i] :
         fmps_lines.append( " BV BV1 " + dc.name + "\n" )

   if ( library.max_dcs_per_pos != 1 ) :
       for i in xrange( library.n_positions ) :
          for ndc_var in ndc_vars[ i ] :
             fmps_lines.append( " BV BV1 " + ndc_var.name + "\n" )
       for sclist in stretch_combos :
          for sc in sclist :
             fmps_lines.append( " BV BV1 " + sc.name + "\n" )
   fmps_lines.append( "ENDATA\n" )

   # #fmps_lines.append( " MARK000 'MARKER' 'INTORG'\n" )
   # 
   # varlist_map = {}
   # recover_var_lines = []
   # 
   # dc = optimize_codons.DegenerateCodon()
   # sum_n_useful_dcs = 0
   # for i in xrange( library.n_positions ) :
   #     count = 0
   #     for j in xrange( 1, library.max_dcs_for_pos[i]+1 ) :
   #        for k in xrange( library.max_per_position_error ) :
   #           if library.divmin_for_error_for_n_dcs[i][j][k] != library.infinity :
   #              count += 1
   #     print i+1, count
   #     sum_n_useful_dcs += count
   # 
   #     vars,ilines = columns_for_position( library, i, dc )
   #     for line in ilines :
   #        fmps_lines.append( line )
   #     varlist_map[ i ] = vars
   # fmps_lines.append( "RHS\n" )
   # fmps_lines.append( " RHS1 libsize " + ( "%6f\n" % math.log( diversity_cap ) ))
   # for i in xrange(library.n_positions) :
   #    fmps_lines.append( " RHS1 onlyone_" + str(i) + " 1\n" )
   # fmps_lines.append( "BOUNDS\n" )
   # varlines = []
   # for i in xrange(library.n_positions) :
   #    ivars = varlist_map[ i ]
   #    for var in ivars :
   #       fmps_lines.append( " BV BV1 " + var[0] + "\n")
   #       varlines.append( var[0] + " " + var[1] + "\n")
   # fmps_lines.append( "ENDATA\n" )
   # 
   open( ilp_input_prefix + ".fmps", "w" ).writelines( fmps_lines )
   varlines = []
   for i in xrange( library.n_positions ) :
      for dc in dcs[ i ] :
         varlines.append( dc.name + " " + dc.codons + " " + dc.aas + "\n" )
   open( ilp_input_prefix + ".dcs", "w" ).writelines( varlines )
