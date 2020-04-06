package SecondaryDTnAnalysisGUI;

import MonteCarloParticleTracer.ParticleType;
import PlottingAPI.Figure;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.font.TextAttribute;
import java.io.File;
import java.text.DateFormat;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.Map;

public class GUI extends JFrame implements WindowListener {

    // *******************
    // Hard Coded Defaults
    // *******************

    // Default number of particles to use per simulation
    private final Double  DEFAULT_NUM_PARTICLES = 5e4;

    // Default number of CPUs to request for each simulation
    private final Integer DEFAULT_NUM_CPUS = Runtime.getRuntime().availableProcessors() - 1;

    // Default exponent to use for temperature profiles T(r) = T0 (1 - r/R)^(power)
    // Uniform model -> 0.0
    private final Double  DEFAULT_PROFILE_EXPONENT = 0.0;

    // Some predefined species for convenience
    private final FuelSpecies[] PREDEFINED_FUEL_SPECIES = new FuelSpecies[] {
            new FuelSpecies("D", 1, 2.014102, 2*2.014102),
            new FuelSpecies("3He", 2, 3.0160293, 3.0160293),
            new FuelSpecies("4He", 2, 4.002602, 4.002602),
            new FuelSpecies("Ar", 18, 39.948, 39.948),
            new FuelSpecies("Kr", 36, 83.798, 83.798),
    };



    // ********************
    // GUI Swing Components
    // ********************

    // Panel that holds all of the Components in this frame
    private JPanel windowPanel = new JPanel(new GridBagLayout());

    // Constraints that we'll edit/use when placing the Components
    private GridBagConstraints constraints = new GridBagConstraints();


    // TextField that holds the primary yield value
    private JFormattedTextField primaryYieldValueField =
            new JFormattedTextField(new DecimalFormat("0.00##E0"));

    // TextField that holds the primary yield uncertainty
    private JFormattedTextField primaryYieldUncertaintyField =
            new JFormattedTextField(new DecimalFormat("0.00##E0"));


    // TextField that holds the secondary yield value
    private JFormattedTextField secondaryYieldValueField =
            new JFormattedTextField(new DecimalFormat("0.00##E0"));

    // TextField that holds the secondary yield uncertainty
    private JFormattedTextField secondaryYieldUncertaintyField =
            new JFormattedTextField(new DecimalFormat("0.00##E0"));


    // TextField that holds the electron temperature value
    private JFormattedTextField electronTemperatureValueField =
            new JFormattedTextField(new DecimalFormat("0.00"));

    // TextField that holds the electron temperature uncertainty
    private JFormattedTextField electronTemperatureUncertaintyField =
            new JFormattedTextField(new DecimalFormat("0.00"));


    // TextField that holds the value of the capsule's outer radius
    private JFormattedTextField outerRadiusValueField =
            new JFormattedTextField(new DecimalFormat("0.00"));

    // TextField that holds the value of the capsule's shell thickness (needed for capsule CR)
    private JFormattedTextField thicknessValueField =
            new JFormattedTextField(new DecimalFormat("0.00"));

    // TextField that holds the value of the capsule's initial fill density TODO User is more likely to know T and P
    private JFormattedTextField initialDensityValueField =
            new JFormattedTextField(new DecimalFormat("0.00"));


    // TextField that holds the number of particles to be used per simulation
    private JFormattedTextField numberParticlesValueField =
            new JFormattedTextField(new DecimalFormat("0.00##E0"));

    // TextField that holds the number of CPUs that will be requested per simulation
    private JFormattedTextField numberCPUsValueField =
            new JFormattedTextField(new DecimalFormat("0"));

    // TextField that holds the exponent to be used in the temperature profiles
    // Uniform Model -> 0.0
    private JFormattedTextField profileExponentValueField =
            new JFormattedTextField(new DecimalFormat("0.00"));


    // Button that lets the user add a species to the plasma model
    private JButton addFuelSpeciesButton = new JButton("Add Species");

    // Button that removes the latest species added to the plasma model
    private JButton removeFuelSpeciesButton = new JButton("Remove Species");

    // Button that starts the analysis
    private JButton startAnalysisButton = new JButton("Start Analysis");


    // Classes that hold all of the TextFields that correspond to individual species
    private ArrayList<FuelSpeciesEntry> fuelSpeciesEntries = new ArrayList<>();


    // Text area that we'll write all logs to to keep track of errors / current status
    private TextArea logConsole = new TextArea();



    /**
     * Main function
     * @param args
     */
    public static void main(String ... args){
        new GUI();
    }


    /**
     * Default constructor called by <code>main()</code>
     */
    private GUI(){
        super("Secondary DTn Analysis GUI");
        setupActionListeners();
        buildMainWindow();
        loadData(Temp_Database.data[0]);
    }


    /**
     * Method that sets ActionListeners to internal methods
     */
    private void setupActionListeners(){

        addFuelSpeciesButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent actionEvent) {
                addFuelSpeciesEntry();
            }
        });

        removeFuelSpeciesButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent actionEvent) {
                removeFuelSpeciesEntry();
            }
        });

        startAnalysisButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent actionEvent) {
                startAnalysis();
            }
        });

    }


    /**
     * Method that constructs this GUI
     */
    private void buildMainWindow(){

        // Temp variables we'll use to keep track of component locations
        int xPos = 0, yPos = 0;


        // (1) First Label
        xPos = 0;
        setConstraints(xPos, yPos, 5, 1);
        setPadding(10, 10, 10, 10);
        windowPanel.add(underlinedJLabel("Measured Values", JLabel.CENTER), constraints);
        yPos++;


        // (2) Primary Yield Line
        xPos = 0;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(5,5,5,5);
        windowPanel.add(new JLabel("Y1n", JLabel.CENTER), constraints);

        xPos++;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(5,5,5,5);
        primaryYieldValueField.setColumns(8);
        primaryYieldValueField.setHorizontalAlignment(JTextField.CENTER);
        primaryYieldValueField.setValue(1e12);
        windowPanel.add(primaryYieldValueField, constraints);

        xPos++;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(10, 10, 10, 10);
        windowPanel.add(new JLabel(" + / - ", JLabel.CENTER), constraints);

        xPos++;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(5,5,5,5);
        primaryYieldUncertaintyField.setColumns(8);
        primaryYieldUncertaintyField.setHorizontalAlignment(JTextField.CENTER);
        primaryYieldUncertaintyField.setValue(1e11);
        windowPanel.add(primaryYieldUncertaintyField, constraints);
        yPos++;


        // (3) Secondary Yield Line
        xPos = 0;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(5,5,5,5);
        windowPanel.add(new JLabel("Y2n", JLabel.CENTER), constraints);

        xPos++;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(5,5,5,5);
        secondaryYieldValueField.setColumns(8);
        secondaryYieldValueField.setHorizontalAlignment(JTextField.CENTER);
        secondaryYieldValueField.setValue(1e10);
        windowPanel.add(secondaryYieldValueField, constraints);

        xPos++;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(5,5,5,5);
        windowPanel.add(new JLabel(" + / - ", JLabel.CENTER), constraints);

        xPos++;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(5,5,5,5);
        secondaryYieldUncertaintyField.setColumns(8);
        secondaryYieldUncertaintyField.setHorizontalAlignment(JTextField.CENTER);
        secondaryYieldUncertaintyField.setValue(1e9);
        windowPanel.add(secondaryYieldUncertaintyField, constraints);
        yPos++;


        // (4) Electron Temperature Line
        xPos = 0;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(5,5,5,5);
        windowPanel.add(new JLabel("Te (keV)", JLabel.CENTER), constraints);

        xPos++;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(5,5,5,5);
        electronTemperatureValueField.setColumns(8);
        electronTemperatureValueField.setHorizontalAlignment(JTextField.CENTER);
        electronTemperatureValueField.setValue(3.0);
        windowPanel.add(electronTemperatureValueField, constraints);

        xPos++;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(5,5,5,5);
        windowPanel.add(new JLabel(" + / - ", JLabel.CENTER), constraints);

        xPos++;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(5,5,5,5);
        electronTemperatureUncertaintyField.setColumns(8);
        electronTemperatureUncertaintyField.setHorizontalAlignment(JTextField.CENTER);
        electronTemperatureUncertaintyField.setValue(0.3);
        windowPanel.add(electronTemperatureUncertaintyField, constraints);
        yPos++;


        // (5) Second Label "Capsule Parameters" & "Model Parameters"
        xPos = 0;
        setConstraints(xPos, yPos, 2, 1);
        setPadding(10, 10, 10, 10);
        windowPanel.add(underlinedJLabel("Capsule Parameters", JLabel.CENTER), constraints);

        xPos = 3;
        setConstraints(xPos, yPos, 2, 1);
        setPadding(10, 10, 10, 10);
        windowPanel.add(underlinedJLabel("Modeling Parameters", JLabel.CENTER), constraints);
        yPos++;


        // (6) Outer Radius Line / Number Particles Line
        xPos = 0;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(5,5,5,5);
        windowPanel.add(new JLabel("Ro (\u03bcm)", JLabel.CENTER), constraints);

        xPos++;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(5,5,5,5);
        outerRadiusValueField.setColumns(8);
        outerRadiusValueField.setHorizontalAlignment(JTextField.CENTER);
        outerRadiusValueField.setValue(1000.0);
        windowPanel.add(outerRadiusValueField, constraints);

        xPos = 3;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(5,5,5,5);
        windowPanel.add(new JLabel("Num Particles", JLabel.CENTER), constraints);

        xPos++;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(5,5,5,5);
        numberParticlesValueField.setColumns(8);
        numberParticlesValueField.setHorizontalAlignment(JTextField.CENTER);
        numberParticlesValueField.setValue(DEFAULT_NUM_PARTICLES);
        windowPanel.add(numberParticlesValueField, constraints);
        yPos++;


        // (7) Thickness / Number of CPUs Line
        xPos = 0;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(5,5,5,5);
        windowPanel.add(new JLabel("t (\u03bcm)", JLabel.CENTER), constraints);

        xPos++;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(5,5,5,5);
        thicknessValueField.setColumns(8);
        thicknessValueField.setHorizontalAlignment(JTextField.CENTER);
        thicknessValueField.setValue(100.0);
        windowPanel.add(thicknessValueField, constraints);

        xPos = 3;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(5,5,5,5);
        windowPanel.add(new JLabel("Num CPUs", JLabel.CENTER), constraints);

        xPos++;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(5,5,5,5);
        numberCPUsValueField.setColumns(8);
        numberCPUsValueField.setHorizontalAlignment(JTextField.CENTER);
        numberCPUsValueField.setValue(DEFAULT_NUM_CPUS);
        windowPanel.add(numberCPUsValueField, constraints);
        yPos++;


        // (8) Initial Pressure Line / Profile Exponent Line
        xPos = 0;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(5,5,5,5);
        windowPanel.add(new JLabel("\u03c1 (mg/cc)", JLabel.CENTER), constraints);

        xPos++;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(5,5,5,5);
        initialDensityValueField.setColumns(8);
        initialDensityValueField.setHorizontalAlignment(JTextField.CENTER);
        initialDensityValueField.setValue(4.0);
        windowPanel.add(initialDensityValueField, constraints);

        xPos = 3;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(5,5,5,5);
        windowPanel.add(new JLabel("Profile Exponent", JLabel.CENTER), constraints);

        xPos++;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(5,5,5,5);
        profileExponentValueField.setColumns(8);
        profileExponentValueField.setHorizontalAlignment(JTextField.CENTER);
        profileExponentValueField.setValue(DEFAULT_PROFILE_EXPONENT);
        windowPanel.add(profileExponentValueField, constraints);
        yPos++;


        // (9) Third Label "Fuel Species"
        xPos = 0;
        setConstraints(xPos, yPos, 5, 1);
        setPadding(10, 10, 10, 10);
        windowPanel.add(underlinedJLabel("Fuel Species", JLabel.CENTER), constraints);
        yPos++;


        // (10) Forth Label "Name Z A Mgas f"
        xPos = 0;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(10, 5, 5, 5);
        windowPanel.add(new JLabel("Name", JLabel.CENTER), constraints);

        xPos++;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(5, 5, 5, 5);
        windowPanel.add(new JLabel("Z", JLabel.CENTER), constraints);

        xPos++;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(5, 5, 5, 5);
        windowPanel.add(new JLabel("A", JLabel.CENTER), constraints);

        xPos++;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(5, 5, 5, 5);
        windowPanel.add(new JLabel("Mgas", JLabel.CENTER), constraints);

        xPos++;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(5, 10, 5, 5);
        windowPanel.add(new JLabel("f", JLabel.CENTER), constraints);
        yPos++;


        // (11) Add / Remove Species Buttons Line
        xPos = 0;
        setConstraints(xPos, yPos, 3, 1);
        setPadding(5,5,5,5);
        windowPanel.add(addFuelSpeciesButton, constraints);

        xPos += 2;
        setConstraints(xPos, yPos, 3, 1);
        setPadding(5,5,5,5);
        windowPanel.add(removeFuelSpeciesButton, constraints);
        yPos++;


        // (12) Start Analysis Button Line
        xPos = 0;
        setConstraints(xPos, yPos, 5, 1);
        setPadding(5,5,5,5);
        startAnalysisButton.setHorizontalAlignment(SwingConstants.CENTER);
        windowPanel.add(startAnalysisButton, constraints);
        yPos++;


        // (13) Line with the log console
        xPos = 0;
        setConstraints(xPos, yPos, 5, 1);
        setPadding(15, 15, 15, 15);
        logConsole.setEditable(false);
        windowPanel.add(logConsole, constraints);


        // Add the panel to this GUI
        this.add(windowPanel);


        // Set up the window
        this.setDefaultCloseOperation(DO_NOTHING_ON_CLOSE);
        this.setLocation(100, 100);         // TODO
        this.addWindowListener(this);
        this.setResizable(false);


        // Show the window
        this.pack();
        this.setVisible(true);


        // Log completion
        writeToLog("Initialization complete!");

    }


    private void loadData(Temp_Database.Data data){
        primaryYieldValueField.setValue(data.Y1n);
        primaryYieldUncertaintyField.setValue(data.Y1n_unc);
        secondaryYieldValueField.setValue(data.Y2n);
        secondaryYieldUncertaintyField.setValue(data.Y2n_unc);
        electronTemperatureValueField.setValue(data.Te);
        electronTemperatureUncertaintyField.setValue(data.Te_unc);

        outerRadiusValueField.setValue(data.Ro);
        thicknessValueField.setValue(data.t);
        initialDensityValueField.setValue(data.rho);
    }


    // ***********************************************************
    // Inherited WindowListener methods that we'll use for cleanup
    // ***********************************************************

    public void windowClosing(WindowEvent e) {
        this.setVisible(false);
        this.dispose();
    }



    // **************************************************
    // Inherited WindowListener methods that we don't use
    // **************************************************

    public void windowOpened     (WindowEvent e) { }
    public void windowClosed     (WindowEvent e) { }
    public void windowIconified  (WindowEvent e) { }
    public void windowDeiconified(WindowEvent e) { }
    public void windowActivated  (WindowEvent e) { }
    public void windowDeactivated(WindowEvent e) { }



    // *********************************
    // Methods called by Action Handlers
    // *********************************

    private void addFuelSpeciesEntry(){

        // Make the entry
        FuelSpeciesEntry entry = new FuelSpeciesEntry();
        fuelSpeciesEntries.add(entry);


        // Temp variables for the position of our components
        int xPos = 0;
        int yPos = 10 + fuelSpeciesEntries.size();


        // Add all of the JComponents associated with this entry
        xPos = 0;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(10, 5, 5, 5);
        windowPanel.add(entry.speciesNameComboBox, constraints);

        xPos++;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(5, 5, 5, 5);
        windowPanel.add(entry.ZTextField, constraints);

        xPos++;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(5, 5, 5, 5);
        windowPanel.add(entry.ATextField, constraints);

        xPos++;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(5, 5, 5, 5);
        windowPanel.add(entry.massTextField, constraints);

        xPos++;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(5, 10, 5, 5);
        windowPanel.add(entry.fractionTextField, constraints);
        yPos++;


        // Move all of the buttons down
        windowPanel.remove(addFuelSpeciesButton);
        windowPanel.remove(removeFuelSpeciesButton);
        windowPanel.remove(startAnalysisButton);


        // Add and remove species buttons
        xPos = 0;
        setConstraints(xPos, yPos, 3, 1);
        setPadding(5,5,5,5);
        windowPanel.add(addFuelSpeciesButton, constraints);

        xPos += 2;
        setConstraints(xPos, yPos, 3, 1);
        setPadding(5,5,5,5);
        windowPanel.add(removeFuelSpeciesButton, constraints);
        yPos++;


        // Start analysis button
        xPos = 0;
        setConstraints(xPos, yPos, 5, 1);
        setPadding(5,5,5,5);
        startAnalysisButton.setHorizontalAlignment(SwingConstants.CENTER);
        windowPanel.add(startAnalysisButton, constraints);
        yPos++;


        // Log Console
        xPos = 0;
        setConstraints(xPos, yPos, 5, 1);
        setPadding(15, 15, 15, 15);
        windowPanel.add(logConsole, constraints);


        // Refresh the GUI
        this.pack();
    }

    private void removeFuelSpeciesEntry(){

        // Sanity check
        if (fuelSpeciesEntries.size() == 0)
            return;


        // Get and remove the latest entry
        FuelSpeciesEntry entry = fuelSpeciesEntries.get(fuelSpeciesEntries.size()-1);
        fuelSpeciesEntries.remove(entry);


        // Remove all of the GUI elements associated with this entry
        windowPanel.remove(entry.speciesNameComboBox);
        windowPanel.remove(entry.ZTextField);
        windowPanel.remove(entry.ATextField);
        windowPanel.remove(entry.massTextField);
        windowPanel.remove(entry.fractionTextField);


        // Refresh the GUI
        this.pack();

    }

    private void startAnalysis(){

        // Init the primitives
        ArrayList<ParticleType> fuelSpecies = new ArrayList<>();
        ArrayList<Double>       numberProportions = new ArrayList<>();

        double Y1n = 0.0, Y1n_Unc = 0.0;
        double Y2n = 0.0, Y2n_Unc = 0.0;
        double Te  = 0.0, Te_Unc  = 0.0;
        double Ro  = 0.0;
        double t   = 0.0;
        double rho = 0.0;

        int numParticles = 0, numCPUs = 0;
        double gamma = 0.0;


        try {

            // Handle the fuel species
            boolean deuteriumPresent = false;
            for (FuelSpeciesEntry entry : fuelSpeciesEntries) {

                ParticleType species = new ParticleType(
                        Double.valueOf(entry.ZTextField.getValue().toString()).intValue(),
                        Double.parseDouble(entry.ATextField.getValue().toString()));
                fuelSpecies.add(species);

                if (species.equals(ParticleType.deuteron))  deuteriumPresent = true;

                numberProportions.add(Double.valueOf(entry.fractionTextField.getValue().toString()));

            }

            if (!deuteriumPresent){
                writeToLog("Deuterium is required to run this simulation!");
                return;
            }

            // Convert all of the other fields to primitives
            Y1n = Double.parseDouble(primaryYieldValueField.getValue().toString());
            Y1n_Unc = Double.parseDouble(primaryYieldUncertaintyField.getValue().toString());

            Y2n = Double.parseDouble(secondaryYieldValueField.getValue().toString());
            Y2n_Unc = Double.parseDouble(secondaryYieldUncertaintyField.getValue().toString());

            Te = Double.parseDouble(electronTemperatureValueField.getValue().toString());
            Te_Unc = Double.parseDouble(electronTemperatureUncertaintyField.getValue().toString());

            Ro = Double.parseDouble(outerRadiusValueField.getValue().toString());
            t = Double.parseDouble(thicknessValueField.getValue().toString());
            rho = Double.parseDouble(initialDensityValueField.getValue().toString());

            numParticles = Double.valueOf(numberParticlesValueField.getValue().toString()).intValue();
            numCPUs = Double.valueOf(numberCPUsValueField.getValue().toString()).intValue();
            gamma = Double.parseDouble(profileExponentValueField.getValue().toString());

        }
        catch (Exception e){
            writeToLog("Fill all fields before starting analysis!");
            return;
        }

        NTOF_Trace specSP = new NTOF_Trace(new File("data/NTOF_Traces/SPEC_SP.dat"));



        // Need to cast stuff to final variables
        final double finalY1n     = Y1n;
        final double finalY1n_Unc = Y1n_Unc;
        final double finalY2n     = Y2n;
        final double finalY2n_Unc = Y2n_Unc;
        final double finalTe      = Te;
        final double finalTe_Unc  = Te_Unc;


        // Build the capsule object
        Capsule capsule = new Capsule(Ro, t, rho, fuelSpecies, numberProportions);

        // Create an analyzer
        Analyzer analyzer = new Analyzer(this, capsule, gamma, numParticles, numCPUs);
        analyzer.getSecondaryYield(new NTOF_Trace[] { specSP });


        // Create a separate thread to run the analysis
        Thread thread = new Thread(new Runnable() {
            @Override
            public void run() {
                analyzer.doAnalysis(finalY1n, finalY1n_Unc, finalY2n, finalY2n_Unc, finalTe, finalTe_Unc);
            }
        });


        // Run the analysis
        thread.start();

    }


    // *****************************
    // Convenience wrapper functions
    // *****************************

    private void setConstraints(int x, int y, int width, int height){
        constraints.gridx = x;
        constraints.gridy = y;
        constraints.gridwidth  = width;
        constraints.gridheight = height;
    }

    private void setPadding(int left, int right, int top, int bottom){
        constraints.insets = new Insets(top, left, bottom, right);
    }

    private JLabel underlinedJLabel(String text, int alignment){

        JLabel label = new JLabel(text, alignment);
        Font font = label.getFont();
        Map attributes = font.getAttributes();
        attributes.put(TextAttribute.UNDERLINE, TextAttribute.UNDERLINE_ON);
        label.setFont(font.deriveFont(attributes));
        return label;

    }

    void writeToLog(String message){

        DateFormat dateFormat = new SimpleDateFormat("HH:mm:ss");
        message = dateFormat.format(new Date()) + " - " + message + "\n";
        logConsole.append(message);

    }

    void writeFormattedToLog(String message, Object ... objects){
        writeToLog(String.format(message, objects));
    }

    void setAnalysisEnabled(boolean enabled){
        startAnalysisButton.setEnabled(enabled);
    }



    // ***************************
    // Private convenience classes
    // ***************************

    private class FuelSpeciesEntry {

        JComboBox speciesNameComboBox;

        JFormattedTextField ZTextField =
                new JFormattedTextField(new DecimalFormat("0"));

        JFormattedTextField ATextField =
                new JFormattedTextField(new DecimalFormat("0.##"));

        JFormattedTextField massTextField =
                new JFormattedTextField(new DecimalFormat("0.##"));

        JFormattedTextField fractionTextField =
                new JFormattedTextField(new DecimalFormat("0.###"));

        public FuelSpeciesEntry() {

            speciesNameComboBox = new JComboBox();
            speciesNameComboBox.setEditable(true);
            ((JLabel) speciesNameComboBox.getRenderer()).setHorizontalAlignment(JLabel.CENTER);

            speciesNameComboBox.setPreferredSize(new Dimension(70, 20));
            speciesNameComboBox.addItem("");
            for (FuelSpecies species : PREDEFINED_FUEL_SPECIES){
                speciesNameComboBox.addItem(species.name);
            }

            speciesNameComboBox.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {

                    // Fill out fields if we recognize the name
                    for (FuelSpecies species : PREDEFINED_FUEL_SPECIES){

                        if (speciesNameComboBox.getSelectedItem().equals(species.name)){
                            ZTextField.setValue(species.Z);
                            ATextField.setValue(species.A);
                            massTextField.setValue(species.gasMass);
                        }

                    }

                    // Refresh the GUI
                    pack();

                }
            });

            ZTextField.setColumns(4);
            ZTextField.setHorizontalAlignment(JTextField.CENTER);

            ATextField.setColumns(4);
            ATextField.setHorizontalAlignment(JTextField.CENTER);

            massTextField.setColumns(4);
            massTextField.setHorizontalAlignment(JTextField.CENTER);

            fractionTextField.setColumns(4);
            fractionTextField.setHorizontalAlignment(JTextField.CENTER);
        }
    }

    private class FuelSpecies {

        private String name;
        private Integer Z;
        private Double A;
        private Double gasMass;

        public FuelSpecies(String name, Integer z, Double a, Double gasMass) {
            this.name = name;
            Z = z;
            A = a;
            this.gasMass = gasMass;
        }
    }


}
