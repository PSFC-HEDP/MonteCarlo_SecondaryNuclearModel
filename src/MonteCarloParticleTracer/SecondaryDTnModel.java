package MonteCarloParticleTracer;

public class SecondaryDTnModel extends Model {


    public SecondaryDTnModel(String name){

        // Call the default constructor
        super(name);

        // Secondary DTn model tracks tritons reacting with a D plasma
        setSourceInformation(NuclearReaction.DD_t, Reactivity.DDp_Reactivity);

        // Secondary DTn model seeks to tally DT neutrons
        addNuclearReaction(NuclearReaction.DT_n);
    }

    public double getYieldRatio(){

        // Grab the DTn tallies
        Tally[] tallies = productParticleEnergyTallyMap.get(NuclearReaction.DT_n);

        // Grab the tally corresponding to neutrons birth
        Tally tally = tallies[0];

        // Return the total of this tally (which happens to be the yield ratio)
        return tally.getTotalWeight();

    }
}
