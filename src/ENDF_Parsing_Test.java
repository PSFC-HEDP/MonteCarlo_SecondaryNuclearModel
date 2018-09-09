import MonteCarloParticleTracer.ENDF_Parser;

import java.io.File;
import java.util.Scanner;

public class ENDF_Parsing_Test {

    private final double NUETRON_MASS_AMU = 1.008664;

    public static void main(String ... args) throws Exception{

        ENDF_Parser parser = new ENDF_Parser(new File("data/ENDFVIII.0_n2H.dat"));
        parser.parseENDF(2);

    }

}
