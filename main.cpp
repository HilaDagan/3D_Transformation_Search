#include "DensityMap.h"
#include "Common.h"
#include "TransformationSearch.h"

// program parameters:
// centfiltmavol011_m.mrc P9Fab.pdb 20
// emd_3061.mrc fragA_41_633_tr.pdb 3.4

// density is according to Chimera:
#define MOL_DENSITY 0.0419
//#define MOL_DENSITY 0.0788

// contains rotations in the quaternion format (4 parameters):
#define ANGLES_FILE "angles_576.qua"
//#define ANGLES_FILE "angels_36864.qua"


/**
 * Find the best transformation between the low resolution map density and the
 * high resolution map of the protein (pdb file), achieved by FFT translations.
 */
int main (int argc, char *argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <mrcFile> <pdbFile> <resolution>" << std::endl;
        return 0;
    }

    // read the map density
    DensityMap mapLo(argv[1]); // the map in low resolution

    // read the molecule
    ChemLib chemLib("chem.lib");
    ChemMolecule molHi;
    Common::readChemMolecule(argv[2], molHi, chemLib);

    // change the centroid of mol high
    Vector3 diff = -(molHi.centroid() - mapLo.getMapCentroid(MOL_DENSITY));
    molHi += diff;

    // Code for debugging and Chimera
//    std::ofstream outFile("molHi_centroid.pdb");
//    outFile << molHi;

    auto resolution = static_cast<float>(atof(argv[3]));

    TransformationSearch transSearch(mapLo, molHi, resolution, mapLo.nx(), mapLo.ny(), mapLo.nz());
    transSearch.getAnglesFromFile(ANGLES_FILE);

    // find the K best translations
    BestK bk = transSearch.startThreaded();

    // print the K best results
    int index = 0;
    for (transformation* trans : bk) {
        std::cerr << "Score: " << trans->score << std::endl;
        std::cerr << "x " << trans->x  << " y " << trans->y << " z " << trans->z
                  << " Angles: " << trans->angles << std::endl;

        ChemMolecule molHiRotatedBack(molHi); // rotate back mol Hi according to the best K
        molHiRotatedBack.rigidTrans(RigidTrans3(trans->angles, Vector3(trans->x, trans->y, trans->z)));
        std::ofstream outFile("out_rotated_"+ std::to_string(index) +".pdb");
        outFile << molHiRotatedBack;
        index++;

        // create map in order to present in Chimera
        DensityMap mapRotated(mapLo);
        mapRotated.simulateMap(molHiRotatedBack, resolution); // create a map density according to the given mol
        mapRotated.writeMap("out_map_rotated"+ std::to_string(index) +".mrc"); // the map after the transformation.
    }
    return 0;
}
