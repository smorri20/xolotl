/**
 * 
 */
package gov.ornl.xolotl.preprocessor;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Properties;

import uk.co.flamingpenguin.jewel.cli.ArgumentValidationException;
import uk.co.flamingpenguin.jewel.cli.CliFactory;
import ncsa.hdf.object.FileFormat;
import ncsa.hdf.object.h5.H5File;

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
	public static void main(String[] args) {

		// Local Declarations
		Arguments myArgs = null;

		// Get command line arguments
		try {
			myArgs = CliFactory.parseArguments(Arguments.class, args);
		} catch (ArgumentValidationException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}

		// Create the Preprocessor - FIXME! Check myArgs != null
		Preprocessor preprocessor = new Preprocessor(myArgs);

		// Generate the network of clusters
		ArrayList<Cluster> clusters = preprocessor.generateNetwork(args);
		
		// Create the HDF5 file
		preprocessor.createHDF5("test.h5");
		
		// Write the header in it
		int[] dim = {8};
		int[] refinement = {0};
		preprocessor.writeHeader("test.h5", dim, refinement);
		
		// Write the network in it
		preprocessor.writeNetwork("test.h5", clusters);

		// Generate the parameters that will be passed to Xolotl
		Properties xolotlParams = preprocessor.generateParameters();
		OutputStream paramsFile = null;

		try {
			// Create the file containing the parameters
			paramsFile = new FileOutputStream("params.txt");

			// Write the parameters to the output file and save 
			// the file to the project root folder
			xolotlParams.store(paramsFile, null);

		} catch (IOException io) {
			io.printStackTrace();
		} finally {
			if (paramsFile != null) {
				try {
					// Close the parameter file
					paramsFile.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}

		return;
	}

}
