import MonteCarloParticleTracer.*;

import java.io.File;
import java.util.ArrayList;
import java.util.Scanner;

public class ModelForNIF3He3He {

    /*
    private static Capsule N170212_003_t1 = new Capsule("D_Nuc_ExPsh_He3_S07a_(t = 2.455 ns)",
            "N170212-003-999", new File("./data/ns1_profile_1d_calb2_t1.dat"), 100);
    private static Capsule N170212_003_t2 = new Capsule("D_Nuc_ExPsh_He3_S07a_(t = 2.555 ns)",
            "N170212-003-999", new File("./data/ns1_profile_1d_calb2_t2.dat"), 100);
    private static Capsule N170212_003_t3 = new Capsule("D_Nuc_ExPsh_He3_S07a_(t = 2.655 ns)",
            "N170212-003-999", new File("./data/ns1_profile_1d_calb2_t3.dat"), 100);
    private static Capsule N170212_003_t4 = new Capsule("D_Nuc_ExPsh_He3_S07a_(t = 2.755 ns)",
            "N170212-003-999", new File("./data/ns1_profile_1d_calb2_t4.dat"), 100);
    private static Capsule N170212_003_t5 = new Capsule("D_Nuc_ExPsh_He3_S07a_(t = 2.855 ns)",
            "N170212-003-999", new File("./data/ns1_profile_1d_calb2_t5.dat"), 100);

    private static Capsule N170212_004_t1 = new Capsule("D_Nuc_ExPsh_He3_S08a_(t = 2.3 ns)",
            "N170212-004-999", new File("./data/ns2_profile_1d_calb2_t1.dat"), 100);
    private static Capsule N170212_004_t2 = new Capsule("D_Nuc_ExPsh_He3_S08a_(t = 2.4 ns)",
            "N170212-004-999", new File("./data/ns2_profile_1d_calb2_t2.dat"), 100);
    private static Capsule N170212_004_t3 = new Capsule("D_Nuc_ExPsh_He3_S08a_(t = 2.5 ns)",
            "N170212-004-999", new File("./data/ns2_profile_1d_calb2_t3.dat"), 100);
    private static Capsule N170212_004_t4 = new Capsule("D_Nuc_ExPsh_He3_S08a_(t = 2.6 ns)",
            "N170212-004-999", new File("./data/ns2_profile_1d_calb2_t4.dat"), 100);
    */

    private static NIF_Shot N170212_003_t1 = new NIF_Shot("ARES_Profiles_2.7ns",
            "N170212-003-999", new File("./data/He3-003_profile_1D_2.7ns.txt"), 120);
    private static NIF_Shot N170212_003_t2 = new NIF_Shot("ARES_Profiles_2.8ns",
            "N170212-003-999", new File("./data/He3-003_profile_1D_2.8ns.txt"), 120);
    private static NIF_Shot N170212_003_t3 = new NIF_Shot("ARES_Profiles_2.9ns",
            "N170212-003-999", new File("./data/He3-003_profile_1D_2.9ns.txt"), 120);
    private static NIF_Shot N170212_003_t4 = new NIF_Shot("ARES_Profiles_3.0ns",
            "N170212-003-999", new File("./data/He3-003_profile_1D_3.0ns.txt"), 120);
    private static NIF_Shot N170212_003_t5 = new NIF_Shot("ARES_Profiles_3.1ns",
            "N170212-003-999", new File("./data/He3-003_profile_1D_3.1ns.txt"), 120);

    private static NIF_Shot N170212_004_t1 = new NIF_Shot("ARES_Profiles_2.5ns",
            "N170212-004-999", new File("./data/He3-004_profile_1D_2.5ns.txt"), 120);
    private static NIF_Shot N170212_004_t2 = new NIF_Shot("ARES_Profiles_2.6ns",
            "N170212-004-999", new File("./data/He3-004_profile_1D_2.6ns.txt"), 120);
    private static NIF_Shot N170212_004_t3 = new NIF_Shot("ARES_Profiles_2.7ns",
            "N170212-004-999", new File("./data/He3-004_profile_1D_2.7ns.txt"), 120);
    private static NIF_Shot N170212_004_t4 = new NIF_Shot("ARES_Profiles_2.8ns",
            "N170212-004-999", new File("./data/He3-004_profile_1D_2.8ns.txt"), 120);

    public static void main(String ... args) throws Exception{

	File cStopPowFile = new File("src/cStopPow/libcStopPow.so");        // Linux
        //File cStopPowFile = new File("src/cStopPow/cStopPow.DLL");       // Windows
        System.load(cStopPowFile.getAbsolutePath());


        NIF_Shot[] shots = {N170212_003_t1, N170212_003_t2, N170212_003_t3, N170212_003_t4, N170212_003_t5,
                N170212_004_t1, N170212_004_t2, N170212_004_t3, N170212_004_t4};

        for (NIF_Shot shot : shots){
            shot.setFill(ParticleType.helium3);
            System.out.println(shot.name);
            runD3HeModel(shot, (int) 1e7, 2);
        }

    }

    public static void runSecondaryDTModel(NIF_Shot shot, Integer numParticles, Integer numCPUs) throws Exception{
        runModel(shot, NuclearReaction.DD_t, Reactivity.DDp_Reactivity, numParticles, numCPUs, NuclearReaction.DT_n);
    }

    public static void runD3HeModel(NIF_Shot shot, Integer numParticles, Integer numCPUs) throws Exception{
        runModel(shot, NuclearReaction.D3He_p, Reactivity.D3Hep_Reactivity, numParticles, numCPUs);
    }

    /*
    public static void run3He3HeModel(Capsule shot, Integer numParticles, Integer numCPUs) throws Exception{

        // Particle distribution

        Distribution energyDistribution  = getWRFDistribution(new File("data/N170212-002-999_10cm_WRF_Spectrum.dat"));
        Distribution burnDistribution = shot.getFuelLayer().getSpatial3He3HepBurnDistribution();

        ParticleDistribution He3He3Distribution = new ParticleDistribution(ParticleType.proton,
                burnDistribution, energyDistribution);


        runModel(shot, He3He3Distribution, numParticles, numCPUs);
    }
    */

    public static void runModel(NIF_Shot shot, NuclearReaction sourceReaction, Reactivity sourceReactivity, Integer numParticles, Integer numCPUs, NuclearReaction ... reactions) throws Exception{

        // Build the model
        String name = String.format("%s_(%s)_%d.output", shot.getName(), shot.getNumber(), System.currentTimeMillis());
        Model model = new Model(name);
        model.addPlasmaLayer(shot.getFuelLayer());
        model.addPlasmaLayer(shot.getShellLayer());
        model.setSourceInformation(sourceReaction, sourceReactivity);


        // Add any nuclear reactions we care about
        for (NuclearReaction reaction : reactions){
            model.addNuclearReaction(reaction);
        }


        // Run the simulation
        model.runSimulation(numParticles, numCPUs);
    }



    static class NIF_Shot {

        private String name;
        private String number;

        private File simulationFile;

        private ParticleType fill;
        private double fuelP0;
        private double totalFuelMass;

        private double[] fuelRadiusNodes;
        private double[] fuelIonTemperatures;
        private double[] fuelEleTemperatures;
        private double[] fuelDensities;

        private double[] shellRadiusNodes;
        private double[] shellIonTemperatures;
        private double[] shellEleTemperatures;
        private double[] shellDensities;

        public NIF_Shot(String name, String number, File simulationFile, int fuelShellInterface){
            this.name = name;
            this.number = number;
            this.simulationFile = simulationFile;

            try {
                Scanner s = new Scanner(simulationFile);
                while (!s.hasNextDouble()) {
                    s.next();
                }

                ArrayList<Double> radiusArray = new ArrayList<>();
                ArrayList<Double> tionArray = new ArrayList<>();
                ArrayList<Double> teArray = new ArrayList<>();
                ArrayList<Double> rhoArray = new ArrayList<>();

                boolean atFuel = true;
                int cellNumber = 0;
                while (s.hasNext()) {
                    cellNumber++;

                    double r = 1e-4 * s.nextDouble();
                    double Tion     = s.nextDouble();
                    double Te       = s.nextDouble();
                    double rho      = s.nextDouble();

                    // This means we hit the boundary
                    if (atFuel && cellNumber > fuelShellInterface) {
                        atFuel = false;

                        fuelRadiusNodes = convertArrayList(radiusArray);
                        fuelDensities = convertArrayList(rhoArray);
                        fuelIonTemperatures = convertArrayList(tionArray);
                        fuelEleTemperatures = convertArrayList(teArray);

                        radiusArray = new ArrayList<>();
                        tionArray = new ArrayList<>();
                        teArray = new ArrayList<>();
                        rhoArray = new ArrayList<>();
                    }

                    // If the first node isn't 0, we'll make a fake zero'th node
                    if (cellNumber == 1 && r != 0){
                        radiusArray.add(0.0);
                        tionArray.add(Tion);
                        teArray.add(Te);
                        rhoArray.add(rho);
                    }

                    // No zero values
                    if (Tion*Te*rho != 0) {
                        radiusArray.add(r);
                        tionArray.add(Tion);
                        teArray.add(Te);
                        rhoArray.add(rho);
                    }
                }

                shellRadiusNodes = convertArrayList(radiusArray);
                shellDensities = convertArrayList(rhoArray);
                shellIonTemperatures = convertArrayList(tionArray);
                shellEleTemperatures = convertArrayList(teArray);

                s.close();
            }
            catch (Exception e){
                e.printStackTrace();
            }

            fuelP0 = fuelRadiusNodes[fuelRadiusNodes.length -1];
            totalFuelMass = getFuelLayer().getTotalMass();

        }

        public void setFill(ParticleType fill) {
            this.fill = fill;
        }

        public void setFuelP0(double fuelP0) {
            this.fuelP0 = fuelP0;
        }

        public void setTotalFuelMass(double totalFuelMass) {
            this.totalFuelMass = totalFuelMass;
        }

        public double getFuelP0() {
            return fuelP0;
        }

        public double getTotalFuelMass() {
            return totalFuelMass;
        }

        public Plasma getFuelLayer(){

            Plasma layer = new Plasma(
                    normalizeNodes(fuelRadiusNodes), fuelIonTemperatures,
                    fuelEleTemperatures, fuelDensities);
            layer.setOuterP0(fuelP0);
            layer.addSpecies(fill, 1.0);

            return layer;
        }

        public Plasma getShellLayer(){

            Plasma layer = new Plasma(
                    normalizeNodes(shellRadiusNodes), shellIonTemperatures,
                    shellEleTemperatures, shellDensities);
            layer.setOuterP0(shellRadiusNodes[shellRadiusNodes.length -1]);
            layer.addSpecies(ParticleType.carbon, 0.5);
            layer.addSpecies(ParticleType.proton, 0.5);

            return layer;
        }

        public String getName() {
            return name;
        }

        public String getNumber() {
            return number;
        }

        public double getCenterTion(){
            return fuelIonTemperatures[0];
        }

    }



    private static double[] convertArrayList(ArrayList<Double> list){
        double[] array = new double[list.size()];

        for (int i = 0; i < list.size(); i++){
            array[i] = list.get(i);
        }

        return array;
    }

    private static double[] normalizeNodes(double[] nodes){
        double[] normNodes = new double[nodes.length];

        double a = nodes[0];
        double b = nodes[nodes.length-1];
        for (int i = 0; i < nodes.length; i++){
            normNodes[i] = (nodes[i] - a) / (b - a);
        }

        return normNodes;
    }

    private static Distribution getWRFDistribution(File file) throws Exception{
        Scanner s = new Scanner(file);

        ArrayList<Double> nodes  = new ArrayList<>();
        ArrayList<Double> values = new ArrayList<>();

        while (s.hasNext()){
            nodes.add(s.nextDouble());
            values.add(s.nextDouble());
        }

        double[] nodeArray  = new double[nodes.size()];
        double[] valueArray = new double[values.size()];

        for (int i = 0; i < nodes.size(); i++){
            nodeArray[i]  =  nodes.get(i);
            valueArray[i] = values.get(i);
        }

        return new Distribution(nodeArray, valueArray);
    }
}
