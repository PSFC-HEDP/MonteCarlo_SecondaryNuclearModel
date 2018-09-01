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

    public void setRadialSourceDistribution(Distribution radialSourceDistribution){
        setSourceInformation(NuclearReaction.DD_t, radialSourceDistribution);
    }

    public double getYieldRatio(){
        return getSecondaryDTNeutronSpectrum().getTotalWeight();
    }

    public Tally getTritonSourceSpectrum(){
        return sourceParticleEnergyTallies[0];
    }

    public Tally getTritionEscapeSpectrum(){
        return sourceParticleEnergyTallies[1];
    }

    public Tally getSecondaryDTNeutronSpectrum(){

        // Grab the DTn tallies
        Tally[] tallies = productParticleEnergyTallyMap.get(NuclearReaction.DT_n);

        // Return the tally corresponding to neutrons birth
        return tallies[0];

    }
}
