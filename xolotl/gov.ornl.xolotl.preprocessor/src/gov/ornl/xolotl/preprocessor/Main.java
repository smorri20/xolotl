/**
 * 
 */
package gov.ornl.xolotl.preprocessor;

import java.util.ArrayList;
import java.util.Properties;

import uk.co.flamingpenguin.jewel.cli.ArgumentValidationException;
import uk.co.flamingpenguin.jewel.cli.CliFactory;

/**
 * This class launches the Xolotl preprocessor.
 * 
 * @author Jay Jay Billings
 * 
 */
public class Main {

	/**
	 * This operation launches the preprocessor and creates the initial
	 * conditions for Xolotl.
	 * 
	 * @param args
	 *            Command line arguments.
	 */
	public static void main(String[] args)
	{

		// Local Declarations
		Arguments myArgs = null;

		// Get command line arguments
		try {	
			myArgs = CliFactory.parseArguments(Arguments.class, args);
		} catch (ArgumentValidationException e1) {
			// TODO Auto-generated catch block
			//e1.printStackTrace();  // Causing error ??
		} 

		// Create the Preprocessor - FIXME! Check myArgs != null
		Preprocessor preprocessor = new Preprocessor(myArgs);

		// Generate the network of clusters
		ArrayList<Cluster> clusters = preprocessor.generateNetwork(args);
		// Dump the clusters to stdout
//		for (Cluster cluster : clusters) {
//			System.out.println(cluster.toString());
//		}
		
		// Create the HDF5 file
		preprocessor.createHDF5("test.h5");
		
		// Write the header in it
		int[] dim = {8};
		int[] refinement = {0};
		preprocessor.writeHeader("test.h5", dim, refinement);
		
		// Write the network in it
		preprocessor.writeNetwork("test.h5", clusters);	

		// Generate the parameters needed to run Xolotl
		Properties xolotlParams = preprocessor.generateParameters();
		// Write the file containing the parameters
		preprocessor.writeParameterFile("params.txt", xolotlParams);
		

		return;
	}				

}
