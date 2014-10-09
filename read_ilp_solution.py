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

if __name__ == "__main__" :
    with blargs.Parser(locals()) as p :
        p.str( "ilp_solution" ).shorthand("s").required()
        p.str( "degenerate_codon_map" ).shorthand( "m" ).required()

    dcm_lines = [ x.strip() for x in open( degenerate_codon_map ).readlines() ]
    ilpsol_lines = [ x.strip() for x in open( ilp_solution ).readlines() ] 

    dcs = { x.split()[0] : ( x.split()[1], x.split()[2] )  for x in dcm_lines }
    found_dc = False
    las_dc = ""
    for line in ilpsol_lines :
        cols = line.split()
        if not found_dc : 
            if len( cols ) < 2 : continue
            if len( cols[1] ) < 3 : continue
            if cols[1][:3] == "DC_" :
                found_dc = True
                last_dc = cols[1]
            else : 
                continue
            if len( cols ) < 6 : continue
            #print line
            assert( cols[2] == "*" )
            if cols[3] == "1" :
                # we found our active DC!
                #print cols
                dc = dcs[ last_dc ]
                print "Active dc:", dc[ 0 ], "aas:", dc[1]
            found_dc = False
            last_dc = ""
        else :
            assert( cols[0] == "*" )
            if cols[1] == "1" :
                # we found our active DC!
                dc = dcs[ last_dc ]
                print "Active dc:", dc[ 0 ], "aas:", dc[1]
            found_dc = False
