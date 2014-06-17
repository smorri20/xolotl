package gov.ornl.xolotl.preprocessor.test;

import static org.junit.Assert.*;

import java.util.Properties;
import java.io.*;

import gov.ornl.xolotl.preprocessor.Preprocessor;

import org.junit.Test;

/**
 * This class is responsible for testing the Preprocessor class
 */
public class PreprocessorTest {

	/**
	 * This operation checks the Helium binding energies.
	 */
	@Test
	public void generateParametersTest() {

		// Local Declarations
		Preprocessor preprocessor = new Preprocessor();
		
		// Generate the default parameters
		Properties defaults = preprocessor.generateParameters();

		// Check that the default material is W
		assertEquals("W", defaults.getProperty("material"));

		// Check that the default startTemp is 1000
		assertEquals("1000", defaults.getProperty("startTemp"));

		// Check that the default tempFile is tempFile
		assertEquals("tempFile", defaults.getProperty("tempFile"));

		// Check that the default heFlux is 2.5e27
		assertEquals("2.5e27", defaults.getProperty("heFlux"));
		
		// Check that the default perfHandler is dummy
		assertEquals("dummy", defaults.getProperty("perfHandler"));
		
		// Check that the default vizHandler is dummy
		assertEquals("dummy", defaults.getProperty("vizHandler"));
		
		// Check that the default perfHandler is dummy
		assertEquals("true", defaults.getProperty("checkpoint"));
		
		// Check that the default perfHandler is dummy
		assertEquals("HDF5file", defaults.getProperty("networkFile"));
		
		// Check that the default perfHandler is dummy
		assertEquals("-ts_final_time 1000 -ts_adapt_dt_max 10 "
				+ "-ts_max_snes_failures 200 -pc_type fieldsplit -pc_fieldsplit_detect_coupling "
				+ "-fieldsplit_0_pc_type redundant -fieldsplit_1_pc_type sor -snes_monitor "
				+ "-ksp_monitor -da_grid_x 10 -ts_max_steps 3 -ts_monitor",
				defaults.getProperty("petscArgs"));

		return;
	}

	/**
	 * This operation makes sure that the parameters can be written
	 * to a file.
	 */
	@Test
	public void testWriteParamsFile() {

		// Local Declarations
		Preprocessor preprocessor = new Preprocessor();
		OutputStream paramsFile = null;
		InputStream inParamsFile = null;
		// Generate the default parameters
		Properties defaults = preprocessor.generateParameters();
		Properties params = new Properties();
		
		try {
			// Create the file containing the parameters
			paramsFile = new FileOutputStream("params.txt");

			// Write the parameters to the file
			defaults.store(paramsFile, null);
			paramsFile.close();

			inParamsFile = new FileInputStream("params.txt");
			// Load the properties file
			params.load(inParamsFile);
			
			// Check that the parameters from the input file are
			// the default parameters
			assertEquals("W", params.getProperty("material"));
			assertEquals("2.5e27", params.getProperty("heFlux"));
			assertEquals("dummy", params.getProperty("perfHandler"));
			assertEquals("1000", params.getProperty("startTemp"));
			
			inParamsFile.close();

		} catch (IOException e) {
			// Complain and fail
			e.printStackTrace();
			fail();
		}

		return;
	}

}
