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

		// Create a file from the uri
		File file = new File("test.hdf5");

		// Retrieve an instance of the HDF5 format
		FileFormat fileFormat = FileFormat
				.getFileFormat(FileFormat.FILE_TYPE_HDF5);

		// Create an H5 file. If it exists already, then delete it.
		try {
			H5File h5File = (H5File) fileFormat.createFile(file.getPath(),
					FileFormat.FILE_CREATE_DELETE);
		} catch (Exception e) {
			// Complain
			e.printStackTrace();
		}

		// Create the Preprocessor - FIXME! Check myArgs != null
		Preprocessor preprocessor = new Preprocessor(myArgs);

		// Generate the network of clusters
		ArrayList<Cluster> clusters = preprocessor.generateNetwork(args);

		// Dump the clusters to stdout
		for (Cluster cluster : clusters) {
			System.out.println(cluster.toString());
		}

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
