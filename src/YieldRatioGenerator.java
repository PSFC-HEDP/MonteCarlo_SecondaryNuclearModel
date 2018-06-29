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

        runModel(Capsule.N180523_001, 0);
        runModel(Capsule.N180523_001, 0.331);
        runModel(Capsule.N180523_001, 0.5);
        runModel(Capsule.N180523_001, 0.75);
        runModel(Capsule.N180523_001, 1.0);


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
            Model model = Model.SecondaryDTnModel(name);
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
            Model model = Model.SecondaryDTnModel(name);
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

    static class Capsule {

        private String shotName;
        private String shotNumber;
        private double burnTemperature;

        private double innerRadius;
        private double outerRadius;

        private double fillDensity;
        private ArrayList<ParticleType> fillSpecies      = new ArrayList<>();
        private ArrayList<Double>       numberProportion = new ArrayList<>();


        private Plasma fuelPlasma;


        public static  Capsule Test_Capsule = new Capsule(
                "",
                "",
                10.0,

                1100,
                100,
                10000,

                new ArrayList<ParticleType>() {{
                    add(ParticleType.deuteron);
                    add(ParticleType.helium3);
                }},

                new ArrayList<Double>(){{
                    add(0.3);
                    add(0.7);
                }}
        );

        public static Capsule N180523_001 = new Capsule(
                "I_Int_Sym_HyC_S02a",
                "N180523-001-999",
                1.85,

                1.1244e+03,
                1.8992e+02,
                6.69,

                new ArrayList<ParticleType>() {{
                    add(ParticleType.deuteron);
                    add(ParticleType.helium3);
                }},

                new ArrayList<Double>(){{
                    add(0.3);
                    add(0.7);
                }}
        );

        public static Capsule N180523_002 = new Capsule(
                "I_Int_Sym_HyC_S01a",
                "N180523-002-999",
                1.65,

                1.1252e+03,
                1.9016e+02,
                6.67,

                new ArrayList<ParticleType>() {{
                    add(ParticleType.deuteron);
                    add(ParticleType.helium3);
                }},

                new ArrayList<Double>(){{
                    add(0.3);
                    add(0.7);
                }}
        );

        public Capsule(String shotName, String shotNumber, double burnTemperature, double outerRadius, double thickness, double fillDensity, ArrayList<ParticleType> fillSpecies, ArrayList<Double> numberProportion) {
            this.shotName           = shotName;
            this.shotNumber         = shotNumber;
            this.burnTemperature    = burnTemperature;

            this.innerRadius        = (outerRadius - thickness) * 1e-4;     // um    -> cm
            this.outerRadius        = outerRadius * 1e-4;                   // um    -> cm
            this.fillDensity        = fillDensity * 1e-3;                   // mg/cc -> g/cc

            this.fillSpecies = fillSpecies;
            this.numberProportion = numberProportion;
        }


        public Plasma getFuelPlasma(double convergence, double gamma){

            // Make the normalized radius nodes
            double[] rNorm = DoubleArray.linspace(0, 1, NUM_SPATIAL_NODES).getValues();


            // We'll assume the density is constant
            double compressedDensity = fillDensity * Math.pow(convergence, 3);
            double[] rho = DoubleArray.linspace(compressedDensity, compressedDensity, NUM_SPATIAL_NODES).getValues();


            // We'll assume a (1 - (r/R)^2)^(gamma) profile
            // Prav is gamma = 0.331
            double[] T = new double[NUM_SPATIAL_NODES];
            for (int i = 0; i < NUM_SPATIAL_NODES; i++){
                T[i] = (1 - Math.pow(rNorm[i], 2));
                T[i] = Math.pow(T[i], gamma);
                T[i] = Math.max(T[i], 0.1);     // We won't let the temperature get too low
            }


            // Calculate the P0
            double P0 = innerRadius / convergence;


            // Make the plasma
            fuelPlasma = new Plasma(rNorm, T, T, rho);
            fuelPlasma.setOuterP0(P0);


            // Add all of the species
            for (int i = 0; i < fillSpecies.size(); i++){
                fuelPlasma.addSpecies(fillSpecies.get(i), numberProportion.get(i));
            }


            // Set the plasma to have the burn temperature
            fuelPlasma.setDDnBurnAveragedIonTemperature(burnTemperature);
            fuelPlasma.setElectronTemperatureFraction(1.0);


            // Return
            return fuelPlasma;
        }

        public double getCapsuleConvergence(double convergence){
            return convergence * (1 + ((outerRadius - innerRadius)/ innerRadius));
        }

        public String getShotName() {
            return shotName;
        }

        public String getShotNumber() {
            return shotNumber;
        }

        public double getInnerRadius() {
            return innerRadius;
        }

        public double getOuterRadius() {
            return outerRadius;
        }

        public double getFillDensity() {
            return fillDensity;
        }
    }



}
