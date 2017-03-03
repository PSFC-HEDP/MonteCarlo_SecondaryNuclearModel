import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;

/**
 * Created by lahmann on 2016-06-15.
 */
public class Particle {

    private ParticleType type;        // This Particle's ID
    private double weight = 1.0;

    private Vector3D position;        // Position vector of this Particle in cm
    private Vector3D velocity;        // Velocity of this Particle in cm/s
    private Vector3D direction;       // Unit vector of the Particle's direction

    private double E;                 // This Particle's energy in MeV




    public Particle(int Z, double mass, Vector3D position, Vector3D direction, double E) {
        this(new ParticleType(Z, mass), position, direction, E);
    }

    public Particle(ParticleType type, Vector3D position, Vector3D direction, double E) {
        this.type = type;
        this.position = position;
        this.direction = direction;
        this.E = E;
    }

    public Particle clone(){
        return new Particle(type, position, direction, E);
    }



    /**
     * Setters
     */

    public void setWeight(double weight) {
        this.weight = weight;
    }

    public void setPosition(Vector3D position) {
        this.position = position;
    }

    public void setDirection(Vector3D direction) {
        this.direction = direction.normalize();
    }



    /**
     * Getters
     */


    public ParticleType getType() {
        return type;
    }

    public double getMass() {
        return type.getMass();
    }

    public int getZ() {
        return type.getZ();
    }

    public double getWeight() {
        return weight;
    }

    public Vector3D getPosition() {
        return position;
    }

    public Vector3D getDirection() {
        return direction;
    }

    public Vector3D getVelocity(){

    }

    public double getE() {
        return E;
    }



    /**
     * Step this Particle some distance through a material with stopping power dEdx
     */

    public Particle step(double distance, double dEdx){
        Vector3D newPosition = this.position.add(distance, direction);
        double newEnergy = this.E + dEdx*distance;

        return new Particle(getType(), newPosition, getDirection(), newEnergy);
    }


    /**
     *
     */

}
