package MonteCarloParticleTracer;

import java.io.File;

class DataFiles {

    // Filenames
    // *********

    private static final String WINDOWS_C_STOP_POW_LIB_FILENAME = "src/cStopPow/cStopPow.DLL";
    private static final String LINUX_C_STOP_POW_LIB_FILENAME   = "src/cStopPow/libcStopPow.so";

    private static final String D3He_XS_FILENAME = "./data/d3He_BoschHale.dat";
    private static final String DT_XS_FILENAME = "./data/dT_BoschHale.dat";

    private static final String DDp_REACTIVITY_FILENAME = "./data/DDp_Reactivities.dat";
    private static final String DDn_REACTIVITY_FILENAME = "./data/DDn_Reactivities.dat";
    private static final String DTn_REACTIVITY_FILENAME = "./data/DTn_Reactivities.dat";
    private static final String D3Hep_REACTIVITY_FILENAME = "./data/D3Hep_Reactivities.dat";
    private static final String He3He3_REACTIVITY_FILENAME = "./data/3He3He_Reactivities.dat";


    private static final String C_NIST_XRAY_MASS_ATTENUATION_FILENAME  = "./data/NIST_xray_C.dat";
    private static final String KAPTON_NIST_XRAY_MASS_ATTENUATION_FILENAME  = "./data/NIST_xray_Kapton.dat";
    private static final String AL_NIST_XRAY_MASS_ATTENUATION_FILENAME = "./data/NIST_xray_Al.dat";
    private static final String V_NIST_XRAY_MASS_ATTENUATION_FILENAME  = "./data/NIST_xray_V.dat";
    private static final String CU_NIST_XRAY_MASS_ATTENUATION_FILENAME = "./data/NIST_xray_Cu.dat";
    private static final String GE_NIST_XRAY_MASS_ATTENUATION_FILENAME = "./data/NIST_xray_Ge.dat";
    private static final String MO_NIST_XRAY_MASS_ATTENUATION_FILENAME = "./data/NIST_xray_Mo.dat";




    // Files
    // *****

    static File windows_StopPow_Lib = new File(WINDOWS_C_STOP_POW_LIB_FILENAME);
    static File linux_StopPow_Lib   = new File(LINUX_C_STOP_POW_LIB_FILENAME);

    static File D3He_XS_File =  new File(D3He_XS_FILENAME);
    static File DT_XS_File   = new File(DT_XS_FILENAME);

    static File DDp_Reactivity_File    = new File(DDp_REACTIVITY_FILENAME);
    static File DDn_Reactivity_File    = new File(DDn_REACTIVITY_FILENAME);
    static File DTn_Reactivity_File    = new File(DTn_REACTIVITY_FILENAME);
    static File D3Hep_Reactivity_File  = new File(D3Hep_REACTIVITY_FILENAME);
    static File He3He3_Reactivity_File = new File(He3He3_REACTIVITY_FILENAME);

    static File C_NIST_MassAttenuationFile  = new File(C_NIST_XRAY_MASS_ATTENUATION_FILENAME);
    static File Kapton_NIST_MassAttenuationFile  = new File(KAPTON_NIST_XRAY_MASS_ATTENUATION_FILENAME);
    static File Al_NIST_MassAttenuationFile = new File(AL_NIST_XRAY_MASS_ATTENUATION_FILENAME);
    static File V_NIST_MassAttenuationFile  = new File(V_NIST_XRAY_MASS_ATTENUATION_FILENAME);
    static File Cu_NIST_MassAttenuationFile = new File(CU_NIST_XRAY_MASS_ATTENUATION_FILENAME);
    static File Ge_NIST_MassAttenuationFile = new File(GE_NIST_XRAY_MASS_ATTENUATION_FILENAME);
    static File Mo_NIST_MassAttenuationFile = new File(MO_NIST_XRAY_MASS_ATTENUATION_FILENAME);

}
