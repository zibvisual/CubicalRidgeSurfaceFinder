// TODO: read mrc files

// read header
// see https://www.ccpem.ac.uk/mrc-format/mrc2014/
// TODO: add explanation of the header

// The unit cell is basically the image. The size of the unit cell (CELLA) is given in angstrom and float (x,y,z).
// MX, MY, MZ are the same as NX,NY,NZ. 

// read data (and convert automatically for now to uint64_t)

// BUT FIRST: debug spatial graph --> debug information is important. Just use any test data for this case --> and save it as numpy file (without space dimensions)