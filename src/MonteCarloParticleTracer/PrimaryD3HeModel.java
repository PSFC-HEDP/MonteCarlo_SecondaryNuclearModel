package MonteCarloParticleTracer;

public class PrimaryD3HeModel extends Model {


    public PrimaryD3HeModel(String name){

        // Call the default constructor
        super(name);

        // Primary D3He-p model just models D3He-protons trying to escape the volume
        setSourceInformation(NuclearReaction.D3He_p, Reactivity.D3Hep_Reactivity);

        // No nuclear reactions are needed for this model

    }

    public Tally getSourceSpectrum(){
        return sourceParticleEnergyTallies[0];
    }

    public Tally getEscapeSpectrum(int layerIndex){
        return sourceParticleEnergyTallies[layerIndex + 1];
    }

}
