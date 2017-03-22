import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.util.FastMath;

/**
 * Created by lahmann on 2016-06-15.
 */
public class Particle {

    private ParticleType type;        // This Particle's ID
    private double weight = 1.0;      // This Particle's statistical weight

    private Vector3D position;        // Position vector of this Particle in cm
    private Vector3D velocity;        // Velocity of this Particle as a fraction of the speed of light
    private Vector3D direction;       // Unit vector of the Particle's direction

    private double energy;             // This Particle's energy in MeV




    public Particle(int Z, double mass, Vector3D position, Vector3D direction, double energy) {
        this(new ParticleType(Z, mass), position, direction, energy);
    }

    public Particle(ParticleType type, Vector3D position, Vector3D direction, double energy) {
        this.type = type;
        this.position = position;
        this.direction = direction;
        this.energy = energy;

        double speed = FastMath.sqrt(2*energy/type.getMass()/Constants.MEV_PER_AMU);        // (v/c)
        this.velocity = direction.scalarMultiply(speed);
    }

    public Particle(ParticleType type, double energy) {
        this.type = type;
        this.position = new Vector3D(0,0,0);
        this.direction = Utils.sampleRandomNormalizedVector();
        this.energy = energy;

        double speed = FastMath.sqrt(2*energy/type.getMass()/Constants.MEV_PER_AMU);        // (v/c)
        this.velocity = direction.scalarMultiply(speed);
    }

    public Particle clone(){
        Particle particle = new Particle(type, position, direction, energy);
        particle.setWeight(weight);
        return particle;
    }



    /**
     * Setters
     */

    public void setWeight(double weight) {
        this.weight = weight;
    }

    public void multiplyWeight(double scalar){
        this.weight *= scalar;
    }

    public void setPosition(Vector3D position) {
        this.position = position;
    }

    public void setDirection(Vector3D direction) {
        this.direction = direction.normalize();
    }

    public void setEnergy(double energy) {
        this.energy = energy;
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
        return velocity;
    }

    public double getEnergy() {
        return energy;
    }



    /**
     * Step this Particle some distance through a material with stopping power dEdx
     */

    public Particle step(double distance, double dEdx){
        Vector3D newPosition = this.position.add(distance, direction);
        double newEnergy = this.energy + dEdx*distance;

        Particle steppedParticle = this.clone();
        steppedParticle.setPosition(newPosition);
        steppedParticle.setEnergy(newEnergy);

        return steppedParticle;
    }

    public String toString(){
        return String.format(
                "Z = %d, M = %.4f\n" +
                "weight = %+.4e\n" +
                "rx = %+.4e cm, ry = %+.4e cm, rz = %+.4e cm\n" +
                "dx = %+.4e   , dy = %+.4e   , dz = %+.4e   \n" +
                "vx = %+.4e   , vy = %+.4e   , vz = %+.4e   \n" +
                "E  = %+.4e MeV\n",

                type.getZ(), type.getMass(),
                weight,
                position.getX(), position.getY(), position.getZ(),
                direction.getX(), direction.getY(), direction.getZ(),
                velocity.getX(), velocity.getY(), velocity.getZ(),
                energy
        );

    }

}
