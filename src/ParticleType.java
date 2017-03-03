/**
 * Created by lahmann on 2017-03-02.
 */
public class ParticleType {

    private int Z;              // Charge in electron charges
    private double mass;        // Mass in amu

    /**
     * Built in pre-defined particles
     */

    public static final ParticleType neutron = new ParticleType(Constants.neutronCharge, Constants.neutronMass_amu);
    public static final ParticleType proton = new ParticleType(Constants.protonCharge, Constants.protonMass_amu);
    public static final ParticleType deuteron = new ParticleType(Constants.deuteriumCharge, Constants.deuteriumMass_amu);
    public static final ParticleType triton = new ParticleType(Constants.tritiumCharge, Constants.tritiumMass_amu);
    public static final ParticleType helium3 = new ParticleType(Constants.helium3Charge, Constants.helium3Mass_amu);
    public static final ParticleType alpha = new ParticleType(Constants.alphaCharge, Constants.alphaMass_amu);


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
}
