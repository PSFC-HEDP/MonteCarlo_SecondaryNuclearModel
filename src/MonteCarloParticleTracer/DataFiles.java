package MonteCarloParticleTracer;

import java.io.File;

class DataFiles {

    // Filenames
    // *********
    private static final String D3He_XS_FILENAME = "./data/d3He_BoschHale.dat";
    private static final String DT_XS_FILENAME = "./data/dT_BoschHale.dat";

    private static final String DDp_REACTIVITY_FILENAME = "./data/DDp_Reactivities.dat";
    private static final String DDn_REACTIVITY_FILENAME = "./data/DDn_Reactivities.dat";
    private static final String DTn_REACTIVITY_FILENAME = "./data/DTn_Reactivities.dat";
    private static final String D3Hep_REACTIVITY_FILENAME = "./data/D3Hep_Reactivities.dat";
    private static final String He3He3_REACTIVITY_FILENAME = "./data/3He3He_Reactivities.dat";



    // Files
    // *****
    static File D3He_XS_File =  new File(D3He_XS_FILENAME);
    static File DT_XS_File   = new File(DT_XS_FILENAME);

    static File DDp_Reactivity_File    = new File(DDp_REACTIVITY_FILENAME);
    static File DDn_Reactivity_File    = new File(DDn_REACTIVITY_FILENAME);
    static File DTn_Reactivity_File    = new File(DTn_REACTIVITY_FILENAME);
    static File D3Hep_Reactivity_File  = new File(D3Hep_REACTIVITY_FILENAME);
    static File He3He3_Reactivity_File = new File(He3He3_REACTIVITY_FILENAME);

}
