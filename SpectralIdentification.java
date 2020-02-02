import java.io.*;
import java.text.DecimalFormat;

public class SpectralIdentification
{
	public static void main(String[] args) throws IOException
	{
		//	Get user input
		System.out.println();
		System.out.println("Element Identification Program Based on Observed Emission Spectrum ");
		System.out.println();
		System.out.print("Please enter an observed spectral line in nanometers (0 to get Table, -1 to quit): ");
		String x1str = IO.readString();

		double x1 = Double.parseDouble(x1str);
		int x2 = (int) x1;
		String[] elements = {"H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca"};
		String spectralDatabaseFilename = "_AtomTrace.asc";
		double[][] lineDB = new double[35][21];
		double dW = 1.0;
		int numElements = 20;										//increase this number to match number of elements in the DB

// Build the default database based on the downloaded emission files.  Other databases can be linked later
// three different sets of Elements data bases are available through https://articles.atomtrace.com/elements-database/
// Elements 1 to 20 from one database were downloaded for this Program
// each file has a complete emission spectrum.
// the program reads each spectrum (file) and finds the wavelengths where there are peaks.
// These are the characterstic emission lines for that element
		for (int i = 1 ; i < 21 ; i++){
			double[][] spectralDatabase = readMatrixFromFile(indexFiler(spectralDatabaseFilename, i));
			for (int j = 0; j < 35; j++) {
				lineDB[j][i-1] = spectralDatabase[j+1][0];
			}
		}

		while (x1 != -1) {
			switch (x2) {

				case 0:

					for ( int j = 0 ; j < numElements ; j++) {
						System.out.print(elements[j] + "\t");
						}
					System.out.println();
					printMatrix(lineDB);
					System.out.print("Please enter an observed spectral line in nanometers (0 to get Table, -1 to quit): ");
					x1str = IO.readString();
					x1 = Double.parseDouble(x1str);
					x2 = (int) x1;
					break;

				default:

					int maxE = 999;
					for (int m = 0; m < 20; m++) {
						for (int n = 0; n < 30; n++) {
							if ((lineDB[n][m] >= (x1 - dW)) && (lineDB[n][m] <= (x1 + dW))) {
								DecimalFormat df = new DecimalFormat("0.0");
								System.out.println();
								System.out.println(elements[m] + " has a spectral line of wavelength: " + df.format(lineDB[n][m]));

								maxE = 0;
								for (int o = 0; o < 30; o++) {
									if (lineDB[o][m] != 0) {
										maxE = o;
									}
								}

								if (maxE != 0) {
									System.out.print("There are " + maxE + " lines for this element: ");
	//								DecimalFormat df = new DecimalFormat("0.0");
									for ( int o = 0; o < maxE; o++) {
										System.out.print(" " + df.format(lineDB[o][m]));
									}
									System.out.println();
								}
							}
						}
					}

					if (maxE == 999){
						System.out.print("No elements found with a spectral line of that wavelength" );
					}


					System.out.println();
					System.out.println();
					System.out.print("Please enter an observed spectral line in nanometers (0 to get Table, -1 to quit): ");
					x1str = IO.readString();
					x1 = Double.parseDouble(x1str);
					x2 = (int) x1;

			}
		}
	}


	public static String indexFiler(String spectralDatabaseFilename, int i) {
		String index = Integer.toString(i);
		String indexFile = index.concat(spectralDatabaseFilename);
		return indexFile;
	}


// This function reads a file ( about 30,000 data points), finds the peak wavlengths and returns these to Main
//The function uses a variation of a function written by Jeff Ames (readMatrixFromFile)
	public static double[][] readMatrixFromFile(String filename) throws IOException
	{
		BufferedReader br;
		String line;
		int numrows, numcols, colsOnLine;
		double[][] matrix;
		String[] lineparts;

		br = new BufferedReader(new FileReader(filename));
		numrows = 0;
		numcols = -1;
		while ((line = br.readLine()) != null)
		{
			numrows++;
			colsOnLine = line.trim().split("\\s+").length;
			if (numcols == -1)
			{
				numcols = colsOnLine;
			}
			if (colsOnLine != numcols)
			{
				System.err.println("Badly formatted matrix file: " + filename);
				return null;
			}
		}
		br.close();

		double iMax = 0;
		double wVal = 0;

// Find largest peak
		br = new BufferedReader(new FileReader(filename));
		for (int row = 0 ; row < numrows ; row++)
		{
			lineparts = br.readLine().trim().split("\\s+");
			double iValue = Double.parseDouble(lineparts[1]);
			if (iValue > iMax) {
				iMax = iValue;
				wVal = Double.parseDouble(lineparts[0]);
			}
		}

		br.close();

// Set threshold based on largest peak
// Find center wavelength of all peaks that exceed that threshold
		double threshold = .05;
		double iThreshold = threshold * iMax;
		double[][] lineSpectrum = new double[40][2];

		int k = 1;
		double a1 = 0;
		double posCross = 0;
		double negCross = 0;
		double peak = 0;

		br = new BufferedReader(new FileReader(filename));
		for (int row = 0 ; row < numrows ; row++)
		{
			lineparts = br.readLine().trim().split("\\s+");

			double iValue = Double.parseDouble(lineparts[1]);
			wVal = Double.parseDouble(lineparts[0]);

			if (a1 < iThreshold && iValue > iThreshold) {
				posCross = wVal;
			}

			if (a1 > iThreshold && iValue < iThreshold) {
				negCross = wVal;
				peak = (posCross + negCross) / 2;				// Peak is in the middle of pos and neg threshold crossing points
				lineSpectrum[k][0] = peak;
				k++;
			}
			a1 = iValue;
		}

		br.close();

		return lineSpectrum;
	}

	public static double integerPart(double d) {
  return (d <= 0) ? Math.ceil(d) : Math.floor(d);
	}

	public static void printMatrix(double[][] matrix)
	{
		DecimalFormat df = new DecimalFormat("0.0");
		int maxrow = 0;
// find number of rows in longest column that has nonzero values

		for (int row = 0 ; row < matrix.length ; row++)
		{
			for (int col = 0 ; col < matrix[0].length - 1  ; col++)
			{
				if (matrix[row][col] != 0) {
					if (row > maxrow) { maxrow = row; }
				}
			}
		}


//Print out table. Don't print zero's

		for ( int row = 0 ; row < maxrow ; row++)
		{
			for ( int col = 0 ; col < matrix[1].length - 1  ; col++)
				{
						if (matrix[row][col] != 0) {
							System.out.print(df.format(matrix[row][col]) + "\t");
						}
					 else{
						 System.out.print("      " + "\t");
					 }
				}
				System.out.println();
	 }

	}

}
