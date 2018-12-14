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

    public static final ParticleType proton = new ParticleType(Constants.HYDROGEN_CHARGE, Constants.PROTON_MASS_AMU);
    public static final ParticleType deuteron = new ParticleType(Constants.HYDROGEN_CHARGE, Constants.DEUTERIUM_MASS_AMU);
    public static final ParticleType triton = new ParticleType(Constants.HYDROGEN_CHARGE, Constants.TRITIUM_MASS_AMU);

    public static final ParticleType helium3 = new ParticleType(Constants.HELIUM_CHARGE, Constants.HELIUM3_MASS_AMU);
    public static final ParticleType alpha = new ParticleType(Constants.HELIUM_CHARGE, Constants.ALPHA_MASS_AMU);

    public static final ParticleType carbon = new ParticleType(Constants.CARBON_CHARGE, Constants.CARBON_MASS_AMU);

    public static final ParticleType tungsten = new ParticleType(Constants.TUNGSTEN_CHARGE, Constants.TUNGSTEN_MASS_AMU);


    public ParticleType(int Z, double mass) {
        this.Z = Z;
        this.mass = mass;
    }

    ParticleType copy(){
        return new ParticleType(Z, mass);
    }

    public boolean equals(ParticleType type){
        return (getZAID() == type.getZAID());
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
        return (getZAID() == type.getZAID());
    }

    // We'll use ZAID as our unique hash identifier
    public int hashCode(){
        return getZAID();
    }

    private int getZAID(){
        return 1000*Z + (int) Math.round(mass);
    }

    public String toString(){
        return String.format("%04d", getZAID());
    }


}
