import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;

/**
 * Created by lahmann on 2016-06-15.
 */
public class Particle {

    private Vector3D position;        // Position vector of this Particle in cm
    private Vector3D direction;       // Unit vector of the Particle's direction
    private double  E;                // This Particle's energy in MeV

    private int A;   // Particle mass number
    private int Z;   // Particle atomic number



    /**
     * Built in pre-defined particles
     */

    public static Particle deuteron(double E){
        return new Particle(E, 1, 2);
    }

    public static Particle triton(double E){
        return new Particle(E, 1, 3);
    }

    public static Particle helium3(double E){
        return new Particle(E, 2, 3);
    }


    /**
     * Constructors
     */

    Particle(double E, int Z, int A) {
        this.E = E;
        this.Z = Z;
        this.A = A;
    }

    public Particle clone(){
        Particle particle = new Particle(this.getE(), this.getZ(), this.getA());
        particle.setPosition(this.getPosition());
        particle.setDirection(this.getDirection());
        return particle;
    }


    /**
     * Setters
     */

    public void setPosition(Vector3D position) {
        this.position = position;
    }

    public void setDirection(Vector3D direction) {
        this.direction = direction.normalize();
    }


    /**
     * Getters
     */

    public Vector3D getPosition() {
        return position;
    }

    public Vector3D getDirection() {
        return direction;
    }

    public double getE() {
        return E;
    }

    public int getA() {
        return A;
    }

    public int getZ() {
        return Z;
    }

    /**
     * Step this Particle some distance through a material with stopping power dEdx
     */

    public Particle step(double distance, double dEdx){
        Vector3D newPosition = this.position.add(distance, direction);
        double newEnergy = this.E + dEdx*distance;

        Particle particle = new Particle(newEnergy, this.getZ(), this.getA());
        particle.setPosition(newPosition);
        particle.setDirection(this.getDirection());

        return particle;
    }


}
