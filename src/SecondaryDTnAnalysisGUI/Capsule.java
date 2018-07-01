package SecondaryDTnAnalysisGUI;

import MonteCarloParticleTracer.DoubleArray;
import MonteCarloParticleTracer.ParticleType;
import MonteCarloParticleTracer.Plasma;

import java.util.ArrayList;

public class Capsule {

    private final int NUM_SPATIAL_NODES = 201;

    private String shotName;
    private String shotNumber;
    private double burnTemperature;

    private double innerRadius;
    private double outerRadius;

    private double fillDensity;
    private ArrayList<ParticleType> fillSpecies      = new ArrayList<>();
    private ArrayList<Double>       numberProportion = new ArrayList<>();



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

    public Capsule(double burnTemperature, double outerRadius, double thickness, double fillDensity, ArrayList<ParticleType> fillSpecies, ArrayList<Double> numberProportion) {
        this("", "", burnTemperature, outerRadius, thickness, fillDensity, fillSpecies, numberProportion);
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
        Plasma fuelPlasma = new Plasma(rNorm, T, T, rho);
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
