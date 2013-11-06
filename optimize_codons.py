import blargs
import genetic_code

class AALibrary :
    def __init__( self ) :
        self.infinity = -1.0;

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
        self.orig_pos = row1[1:]
        self.max_obs = 0
        for i in xrange(20) :
            line = lines[ i + 1 ]
            vals = line.split(",")[1:]
            for j in xrange(len(vals)):
                self.aa_counts[j][i] = int(vals[j])
                if self.aa_counts[j][i] > self.max_obs :
                    self.max_obs = self.aa_counts[j][i]

    def compute_smallest_diversity_for_all_errors( self ) :
        self.divmin_for_error = [ [] ] * self.n_positions
        for i in xrange( self.n_positions ) : self.divmin_for_error[i] = [ self.infinity ] * self.max_obs
        

if __name__ == "__main__" :
    with blargs.Parser( locals() ) as p :
        p.str( "input_csv" ).required()
        p.float( "diversity_cap" ).required()

    library = AALibrary()
    library.load_library( input_csv )
    print library.max_obs
    print library.aa_counts
