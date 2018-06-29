package AnalyzerGUI;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.font.TextAttribute;
import java.text.DecimalFormat;
import java.text.Format;
import java.util.ArrayList;
import java.util.Map;

public class GUI extends JFrame implements WindowListener {


    private JPanel windowPanel = new JPanel(new GridBagLayout());

    private GridBagConstraints constraints = new GridBagConstraints();


    private JFormattedTextField primaryYieldValueField =
            new JFormattedTextField(new DecimalFormat("0.00##E0"));
    private JFormattedTextField primaryYieldUncertaintyField =
            new JFormattedTextField(new DecimalFormat("0.00##E0"));


    private JFormattedTextField secondaryYieldValueField =
            new JFormattedTextField(new DecimalFormat("0.00##E0"));
    private JFormattedTextField secondaryYieldUncertaintyField =
            new JFormattedTextField(new DecimalFormat("0.00##E0"));


    private JFormattedTextField electronTemperatureValueField =
            new JFormattedTextField(new DecimalFormat("0.00"));
    private JFormattedTextField electronTemperatureUncertaintyField =
            new JFormattedTextField(new DecimalFormat("0.00"));


    private JFormattedTextField innerRadiusValueField =
            new JFormattedTextField(new DecimalFormat("0.00"));

    private JFormattedTextField outerRadiusValueField =
            new JFormattedTextField(new DecimalFormat("0.00"));

    private JFormattedTextField initialPressureValueField =
            new JFormattedTextField(new DecimalFormat("0.00"));


    private JButton addFuelSpeciesButton = new JButton("Add Species");


    private ArrayList<FuelSpeciesEntry> fuelSpeciesEntries = new ArrayList<>();



    public static void main(String ... args){
        new GUI();
    }


    private GUI(){
        super("Secondary Analysis GUI");
        initialize();
        mainLoop();
    }


    private void initialize(){

        setupActionListeners();
        buildMainWindow();

    }


    private void setupActionListeners(){

        addFuelSpeciesButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent actionEvent) {
                addFuelSpeciesEntry();
            }
        });

    }


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
        windowPanel.add(primaryYieldValueField, constraints);

        xPos++;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(10, 10, 10, 10);
        windowPanel.add(new JLabel(" +/- ", JLabel.CENTER), constraints);

        xPos++;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(5,5,5,5);
        primaryYieldUncertaintyField.setColumns(8);
        primaryYieldUncertaintyField.setHorizontalAlignment(JTextField.CENTER);
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
        windowPanel.add(secondaryYieldValueField, constraints);

        xPos++;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(5,5,5,5);
        windowPanel.add(new JLabel(" +/- ", JLabel.CENTER), constraints);

        xPos++;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(5,5,5,5);
        secondaryYieldUncertaintyField.setColumns(8);
        secondaryYieldUncertaintyField.setHorizontalAlignment(JTextField.CENTER);
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
        windowPanel.add(electronTemperatureValueField, constraints);

        xPos++;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(5,5,5,5);
        windowPanel.add(new JLabel(" +/- ", JLabel.CENTER), constraints);

        xPos++;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(5,5,5,5);
        electronTemperatureUncertaintyField.setColumns(8);
        electronTemperatureUncertaintyField.setHorizontalAlignment(JTextField.CENTER);
        windowPanel.add(electronTemperatureUncertaintyField, constraints);
        yPos++;


        // (5) Second Label "Capsule Parameters"
        xPos = 0;
        setConstraints(xPos, yPos, 5, 1);
        setPadding(10, 10, 10, 10);
        windowPanel.add(underlinedJLabel("Capsule Parameters", JLabel.CENTER), constraints);
        yPos++;


        // (6) Inner Radius Line
        xPos = 0;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(5,5,5,5);
        windowPanel.add(new JLabel("Ri (\u03bcm)", JLabel.CENTER), constraints);

        xPos++;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(5,5,5,5);
        innerRadiusValueField.setColumns(8);
        innerRadiusValueField.setHorizontalAlignment(JTextField.CENTER);
        windowPanel.add(innerRadiusValueField, constraints);
        yPos++;


        // (7) Outer Radius Line
        xPos = 0;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(5,5,5,5);
        windowPanel.add(new JLabel("Ro (\u03bcm)", JLabel.CENTER), constraints);

        xPos++;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(5,5,5,5);
        outerRadiusValueField.setColumns(8);
        outerRadiusValueField.setHorizontalAlignment(JTextField.CENTER);
        windowPanel.add(outerRadiusValueField, constraints);
        yPos++;


        // (8) Initial Pressure Line
        xPos = 0;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(5,5,5,5);
        windowPanel.add(new JLabel("\u03c1 (mg/cc)", JLabel.CENTER), constraints);

        xPos++;
        setConstraints(xPos, yPos, 1, 1);
        setPadding(5,5,5,5);
        initialPressureValueField.setColumns(8);
        initialPressureValueField.setHorizontalAlignment(JTextField.CENTER);
        windowPanel.add(initialPressureValueField, constraints);
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


        // (11) Add Species Button Line
        xPos = 0;
        setConstraints(xPos, yPos, 5, 1);
        setPadding(10, 10, 10, 10);
        windowPanel.add(addFuelSpeciesButton, constraints);
        yPos++;

        // Add the panel to this GUI
        this.add(windowPanel);


        // Set up the window
        this.setDefaultCloseOperation(DO_NOTHING_ON_CLOSE);
        this.setLocation(600, 600);         // TODO
        this.addWindowListener(this);
        this.setResizable(false);


        // Show the window
        this.pack();
        this.setVisible(true);

    }

    private void mainLoop(){

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


        // Move the add species button down
        windowPanel.remove(addFuelSpeciesButton);

        xPos = 0;
        setConstraints(xPos, yPos, 5, 1);
        setPadding(10, 10, 10, 10);
        windowPanel.add(addFuelSpeciesButton, constraints);


        // Refresh the GUI
        this.pack();
    }

    private JLabel underlinedJLabel(String text, int alignment){

        JLabel label = new JLabel(text, alignment);
        Font font = label.getFont();
        Map attributes = font.getAttributes();
        attributes.put(TextAttribute.UNDERLINE, TextAttribute.UNDERLINE_ON);
        label.setFont(font.deriveFont(attributes));
        return label;

    }

    private class FuelSpeciesEntry {

        JComboBox speciesNameComboBox = new JComboBox();

        JFormattedTextField ZTextField =
                new JFormattedTextField(new DecimalFormat("0"));

        JFormattedTextField ATextField =
                new JFormattedTextField(new DecimalFormat("0"));

        JFormattedTextField massTextField =
                new JFormattedTextField(new DecimalFormat("0"));

        JFormattedTextField fractionTextField =
                new JFormattedTextField(new DecimalFormat("0.###"));

        public FuelSpeciesEntry() {
            speciesNameComboBox.setEditable(true);

            ZTextField.setColumns(4);
            ATextField.setColumns(4);
            massTextField.setColumns(4);
            fractionTextField.setColumns(4);
        }
    }


}
