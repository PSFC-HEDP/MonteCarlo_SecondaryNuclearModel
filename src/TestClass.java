/**
 * Created by lahmann on 2016-06-15.
 * TODO: Test difference between hotspot and uniform and formal profiles
 * TODO: Add spatial distributions based on T and \rho
 * TODO: Add simple ionization in Plasma for getting ne
 * TODO: Move some of the dEdx calculations to Plasma
 * TODO: Add getRhoR to Plasma class (for comparing models)
 * TODO: Add getMeanTe to Plasma class (for comparing models)
 * TODO: Add plotting capabilities (MATLAB)
 */
public class TestClass {

    public static void main(String ... args) throws Exception {

        System.load("C:\\Dropbox (MIT)\\Research Work\\SecondaryAnalysis_Java\\src\\cStopPow\\cStopPow.DLL");

        System.out.print("Generating plasma ... ");
        Plasma plasma = Plasma.uniformPlasma(50.0, 5.0, 3.0);
        plasma.addDeuteriumSpecies(1.0);
        System.out.println("Done!");

        System.out.print("Generating triton distribution ... ");
        ParticleDistribution tritions =
                ParticleDistribution.ThermalDDtDistribution(plasma.getSpatialDDpBurnDistribution(), 1e-3);
        System.out.println("Done!");

        System.out.print("Building model ... ");
        Model model = new Model(tritions, CrossSection.dt(), plasma);
        System.out.println("Done!");

        System.out.print("Running model ... \n");
        System.out.println("Yield ratio = " + model.getYieldRatio((int) 1e4));
        System.out.println("Done!");


    }



}
