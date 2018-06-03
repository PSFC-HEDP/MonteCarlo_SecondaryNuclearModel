package MonteCarloParticleTracer;

import java.io.File;
import java.io.FileWriter;

public class OutputFile extends File {

    private final String SEPERATOR = "********************************************************************************";


    public OutputFile(String s) {
        super(s);
    }

    public void addHeader(String header) throws Exception{

        String centeredHeader = "* ";

        int numSpaces = Math.floorDiv(76 - header.length(), 2);
        for (int i = 0; i < numSpaces; i++){
            centeredHeader += " ";
        }
        centeredHeader += header;

        while (centeredHeader.length() < 79){
            centeredHeader += " ";
        }
        centeredHeader += "*";


        addString(SEPERATOR);
        addString(centeredHeader);
        addString(SEPERATOR);
        addString("");
    }

    public void addSubheader(String subHeader) throws Exception{
        addString("--- " + subHeader + " ---");
        addString("");
    }

    public void addString(String string) throws Exception{
        FileWriter w = new FileWriter(this.getAbsoluteFile(), true);
        w.write(string + "\n");
        w.close();
    }
}
