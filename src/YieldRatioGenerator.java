import SecondaryDTnAnalysisGUI.Capsule;
import MonteCarloParticleTracer.*;

import java.io.File;
import java.util.ArrayList;
import java.util.Scanner;


public class YieldRatioGenerator {

    private static final int NUM_SPATIAL_NODES  = (int) 201;
    private static final int NUM_CPUS           = (int) 47;
    private static final int NUM_PARTICLES      = (int) (NUM_CPUS * 1e3);
    //private static final int NUM_CPUS           = Runtime.getRuntime().availableProcessors() - 1;

    public static void main(String ... args) throws Exception{

        /*
        runModel(Capsule.N180523_001, 0);
        runModel(Capsule.N180523_001, 0.331);
        runModel(Capsule.N180523_001, 0.5);
        runModel(Capsule.N180523_001, 0.75);
        runModel(Capsule.N180523_001, 1.0);
        */


        //uniformModelBenchmark(10.0, 10.0);


    }

    public static void uniformModelBenchmark(double T, double rho)  throws Exception{

        double[] rhoRs = DoubleArray.logspace(-2, 1, 201).getValues();       // g/cc

        for (double rhoR : rhoRs){

            // Make the model name
            String name = String.format("Uniform_Model_%d.output", System.currentTimeMillis());


            // Build the uniform fuel plasma
            double[] rs   = DoubleArray.linspace(0, 1, 201).getValues();
            double[] Ts   = DoubleArray.linspace(T, T, 201).getValues();
            double[] rhos = DoubleArray.linspace(rho, rho, 201).getValues();
            Plasma fuelPlasma = new Plasma(rs, Ts, Ts, rhos);

            double P0 = rhoR / rho;      // cm
            fuelPlasma.setOuterP0(P0);

            fuelPlasma.addSpecies(ParticleType.deuteron, 0.3);
            fuelPlasma.addSpecies(ParticleType.helium3, 0.7);


            // Build the model
            SecondaryDTnModel model = new SecondaryDTnModel(name);
            model.addPlasmaLayer(fuelPlasma);


            // Run the simulation
            model.runSimulation(NUM_PARTICLES, NUM_CPUS);


            // Get the output file
            File output = model.getOutputFile();


            // Really shitty parser
            Scanner s = new Scanner(output);
            ArrayList<Double> totals = new ArrayList<>();
            while (s.hasNext()){
                String line = s.nextLine();

                if (line.contains("total")){
                    String[] temp = line.split(",");
                    totals.add(Double.valueOf(temp[1]));
                }
            }
            s.close();

            double ratio = totals.get(7);
            System.out.printf("%.4e %.4e\n", rhoR, ratio);

        }

    }

    public static void runModel(Capsule capsule, double gamma) throws Exception{

        double[] CRs = DoubleArray.linspace(5.0, 30.0, 101).getValues();

        for (double CR : CRs) {

            // Build the model name
            String name = String.format("%s_(%s)_%d.output",
                    capsule.getShotName(),
                    capsule.getShotNumber(),
                    System.currentTimeMillis());


            // Build the fuel plasma
            Plasma fuelPlasma = capsule.getFuelPlasma(CR, gamma);


            // Build the model
            SecondaryDTnModel model =new SecondaryDTnModel(name);
            model.addPlasmaLayer(fuelPlasma);


            // Run the simulation
            model.runSimulation(NUM_PARTICLES, NUM_CPUS);


            // Get the output file
            File output = model.getOutputFile();


            // Really shitty parser
            Scanner s = new Scanner(output);
            ArrayList<Double> totals = new ArrayList<>();
            while (s.hasNext()){
                String line = s.nextLine();

                if (line.contains("total")){
                    String[] temp = line.split(",");
                    totals.add(Double.valueOf(temp[1]));
                }
            }
            s.close();

            double ratio = totals.get(7);
            System.out.printf("%.2f %.2f %.2f %.4e %.4e\n",
                    CR,
                    capsule.getCapsuleConvergence(CR),
                    10000*capsule.getInnerRadius() / CR,
                    1000*fuelPlasma.getArealDensity(0.0, 0.0),
                    ratio);
        }

    }





}
