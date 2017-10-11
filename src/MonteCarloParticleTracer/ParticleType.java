package MonteCarloParticleTracer;

/**
 * Created by lahmann on 2017-03-02.
 */
public class ParticleType {

    private int Z;              // Charge in electron charges
    private double mass;        // Mass in amu

    /**
     * Built in pre-defined particles
     */

    public static final ParticleType neutron = new ParticleType(Constants.NEUTRON_CHARGE, Constants.NEUTRON_MASS_AMU);
    public static final ParticleType proton = new ParticleType(Constants.PROTON_CHARGE, Constants.PROTON_MASS_AMU);
    public static final ParticleType deuteron = new ParticleType(Constants.DEUTERIUM_CHARGE, Constants.DEUTERIUM_MASS_AMU);
    public static final ParticleType triton = new ParticleType(Constants.TRITIUM_CHARGE, Constants.TRITIUM_MASS_AMU);
    public static final ParticleType helium3 = new ParticleType(Constants.HELIUM3_CHARGE, Constants.HELIUM3_MASS_AMU);
    public static final ParticleType alpha = new ParticleType(Constants.ALPHA_CHARGE, Constants.ALPHA_MASS_AMU);
    public static final ParticleType carbon = new ParticleType(Constants.CARBON_CHARGE, Constants.CARBON_MASS_AMU);


    public ParticleType(int Z, double mass) {
        this.Z = Z;
        this.mass = mass;
    }

    public boolean equals(ParticleType type){
        return (Z == type.getZ() && mass == type.getMass());
    }

    public int getZ() {
        return Z;
    }

    public double getMass() {
        return mass;
    }

    public boolean equals(Object o){
        if (o == this) return true;
        if (!(o instanceof ParticleType)) return false;

        ParticleType type = (ParticleType) o;
        return Z == type.getZ() && mass == type.getMass();
    }

    // We'll use ZAID as our unique hash identifier
    public int hashCode(){
        return 1000*Z + (int) Math.round(mass);
    }

    public String toString(){
        int ZAID = 1000*Z + (int) Math.round(mass);
        return Integer.toString(ZAID);
    }


}