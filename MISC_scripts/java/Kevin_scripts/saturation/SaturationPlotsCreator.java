import com.lowagie.text.*;
import com.lowagie.text.pdf.*;

import org.jfree.chart.*;
import org.jfree.chart.plot.*;
import org.jfree.data.statistics.*;
import org.jfree.ui.*;

import java.awt.*;
import java.awt.geom.*;

import java.io.*;
import java.util.*;

/**
 *	<pre>
 *	This class produces saturation plots from a set of files, each
 *	containing a list of genomic regions. Each genomic region is specified
 *	using the following format:
 *
 *	<ID><tab><start><tab><end>
 *
 *	where
 *	<ID>	is the identifier of the region-at-large, such as the chromosome
 *	<start>	is the starting position of the region
 *	<end>	is the ending position of the region (this position is inside
 *		the region)
 *
 *	The y-axis could be the absolute number of nucleotides, or a fraction of
 *	an input total number of nucleotides, such as the total number of
 *	nucleotides of the coding transcripts in the example. To use the absolute
 *	number, input the total as 0.
 *
 *	If the number of input files is no more than 31, the program can compute
 *	the coverage from all combinations of the input datasets. If the number
 *	of input files is more than 31, or if the number of combinations is more
 *	than a specified threshold, a random sample of the combinations will be
 *	considered.
 *	</pre>
 *
 *	@author		Kevin Yuk-Lap Yip, Gerstein Lab, Yale University
 *			(project with Zhi John Lu)
 *	@version	1.2 (March 5, 2010)
 *
 *	<pre>
 *	Change History:
 *	1.2	- Added option for printing the count/fraction statistics
 *	1.1	- Added option for sampling the combinations
 *	1.0	- Initial version
 *	</pre>
 */

public class SaturationPlotsCreator
{
	/*----------------------------------------------------------------------
	 Constants
	----------------------------------------------------------------------*/

	/**
	 *	Default maximum number of combinations.
	 */
	protected static final int DEF_MAX_COM = 1024;


	/*----------------------------------------------------------------------
	 Public methods
	----------------------------------------------------------------------*/

	public static void main(String[] argv) throws Exception
	{
		// Check the number of arguments, and print usage if necessary
		if (argv.length < 5)
		{
			System.err.println("Usage: java -classpath <classpath> SaturationPlotsCreator <Output graph file (.pdf)> <Output count/fraction statistics file (.txt)> <Total number of nucleotides, 0 for plotting absolute numbers> <Maximum number of combinations, 0 to consider all combinations> <Input file 1> <Input file 2> ...");
			System.exit(-1);
		}

		// Register input parameters
		int argc = 0;
		File outGraphFile = new File(argv[argc++]);
		File outCountFile = new File(argv[argc++]);
		double totNtCount = Double.parseDouble(argv[argc++]);	// Use double to make divisions easy
		int maxCom = Integer.parseInt(argv[argc++]);

		int inFileCount = argv.length - argc;
		File[] inFiles = new File[inFileCount];
		for (int i=0; i<inFileCount; i++)
			inFiles[i] = new File(argv[argc++]);

		// Scan the input files once to get all end points
		Map<String, Set<Integer>> endPoints = new HashMap<String, Set<Integer>>();
		for (int i=0; i<inFileCount; i++)
		{
			BufferedReader br = new BufferedReader(new FileReader(inFiles[i]));
			String line = br.readLine();
			while (line != null && line.length() != 0)
			{
				StringTokenizer st = new StringTokenizer(line, "\t", false);
				String id = st.nextToken();
				Set<Integer> idEndPoints = endPoints.get(id);
				if (idEndPoints == null)
				{
					idEndPoints = new HashSet<Integer>();
					endPoints.put(id, idEndPoints);
				}
				idEndPoints.add(new Integer(st.nextToken()));				// Start
				idEndPoints.add(new Integer(Integer.parseInt(st.nextToken()) + 1));	// End + 1
				line = br.readLine();
			}
			br.close();
		}

		// Store the end points in arrays for more efficient access
		Map<String, int[]> fastEndPoints = new TreeMap<String, int[]>();
		for (String id : endPoints.keySet())
		{
			Set<Integer> idEndPoints = endPoints.get(id);
			Iterator<Integer> iter = idEndPoints.iterator();
			int[] fastIdEndPoints = new int[idEndPoints.size()];
			for (int i=0; i<idEndPoints.size(); i++)
				fastIdEndPoints[i] = iter.next().intValue();
			Arrays.sort(fastIdEndPoints);
			fastEndPoints.put(id, fastIdEndPoints);
		}

		// Start counting
		java.util.List<java.util.List<Double>> ntDists = new ArrayList<java.util.List<Double>>(inFileCount);

		// Case 1: consider all combinations
		if (inFileCount <= 31 && (maxCom == 0 || Math.pow(2, inFileCount) <= maxCom))
		{
			// Phase 1: for each disjoint region, record which features have it
			int[] featureSigs = new int[inFileCount];
			featureSigs[0] = 1;
			for (int i=1; i<inFileCount; i++)
				featureSigs[i] = featureSigs[i-1] << 1;

			Map<String, int[]> coverMap = new HashMap<String, int[]>();
			for (String id : fastEndPoints.keySet())
			{
				int[] coveredFeatures = new int[fastEndPoints.get(id).length - 1];
				coverMap.put(id, coveredFeatures);
			}

			for (int i=0; i<inFileCount; i++)
			{
				System.out.println("Performing counting for " + inFiles[i].getName());
				BufferedReader br = new BufferedReader(new FileReader(inFiles[i]));
				String line = br.readLine();
				while (line != null && line.length() != 0)
				{
					// Get the region
					StringTokenizer st = new StringTokenizer(line, "\t", false);
					String id = st.nextToken();
					int[] fastIdEndPoints = fastEndPoints.get(id);
					int[] coveredFeatures = coverMap.get(id);
					int start = Integer.parseInt(st.nextToken());
					int end = Integer.parseInt(st.nextToken());

					// Do the countings
					int endPointIndex = Arrays.binarySearch(fastIdEndPoints, start);
					while (fastIdEndPoints[endPointIndex] != end + 1)
					{
						if ((coveredFeatures[endPointIndex] & featureSigs[i]) == 0)
							coveredFeatures[endPointIndex] += featureSigs[i];
						endPointIndex ++;
					}

					// Read the next line
					line = br.readLine();
				}
				br.close();
			}

			// Phase 2: for each disjoint region, increment the counters accordingly
			System.out.println("Aggregating the counts.");
			int countNum = (int)Math.pow(2, inFileCount);
			int[] ntCounts = new int[countNum];

			for (String id : fastEndPoints.keySet())
			{
				int[] fastIdEndPoints = fastEndPoints.get(id);
				int[] coveredFeatures = coverMap.get(id);

				for (int i=0; i<fastIdEndPoints.length-1; i++)
				{
					int start = fastIdEndPoints[i];
					int endExc = fastIdEndPoints[i+1];
					int lengthInc = endExc - start;

					for (int j=0; j<countNum; j++)
						if ((j & coveredFeatures[i]) != 0)

							ntCounts[j] += lengthInc;
				}
			}

			// Determine the number of features each counter is related to
			int[] relatedFeatureCounts = new int[countNum];
			for (int i=0; i<inFileCount; i++)
			{
				int ptr1 = 0;
				int ptr2 = 0;
				for (int j=0; j<Math.pow(2, inFileCount-i-1); j++)
				{
					for (int k=0; k<Math.pow(2, i); k++)
						ptr1 ++;
					for (int k=0; k<Math.pow(2, i); k++)
					{
						relatedFeatureCounts[ptr1] ++;
						ptr1 ++;
					}
				}
			}

			// Store counts in lists
			for (int i=0; i<inFileCount; i++)
			{
				int pointCount = choose(inFileCount, i+1);
				ntDists.add(new ArrayList<Double>(pointCount));
			}

			if (totNtCount == 0)
				for (int i=1; i<countNum; i++)
					ntDists.get(relatedFeatureCounts[i] - 1).add(new Double(ntCounts[i]));
			else
				for (int i=1; i<countNum; i++)
					ntDists.get(relatedFeatureCounts[i] - 1).add(new Double(ntCounts[i]/ totNtCount));
		}

		// Case 2: consider a sample of the combinations
		else
		{
			// Use the default maximum number of combinations if not
			// specified
			if (maxCom == 0)
				maxCom = DEF_MAX_COM;
			Random rand = new Random(0);

			// Phase 1: for each disjoint region, record which features have it
			Map<String, BitSet[]> coverMap = new HashMap<String, BitSet[]>();
			for (String id : fastEndPoints.keySet())
			{
				BitSet[] coveredFeatures = new BitSet[fastEndPoints.get(id).length - 1];
				for (int i=0; i<fastEndPoints.get(id).length - 1; i++)
					coveredFeatures[i] = new BitSet(inFileCount);
				coverMap.put(id, coveredFeatures);
			}

			for (int i=0; i<inFileCount; i++)
			{
				System.out.println("Performing counting for " + inFiles[i].getName());
				BufferedReader br = new BufferedReader(new FileReader(inFiles[i]));
				String line = br.readLine();
				while (line != null && line.length() != 0)
				{
					// Get the region
					StringTokenizer st = new StringTokenizer(line, "\t", false);
					String id = st.nextToken();
					int[] fastIdEndPoints = fastEndPoints.get(id);
					BitSet[] coveredFeatures = coverMap.get(id);
					int start = Integer.parseInt(st.nextToken());
					int end = Integer.parseInt(st.nextToken());

					// Do the countings
					int endPointIndex = Arrays.binarySearch(fastIdEndPoints, start);
					while (fastIdEndPoints[endPointIndex] != end + 1)
					{
						coveredFeatures[endPointIndex].set(i);
						endPointIndex ++;
					}

					// Read the next line
					line = br.readLine();
				}
				br.close();
			}

			// Phase 2: for each disjoint region, increment the counters accordingly
			System.out.println("Aggregating the counts.");
			for (int fileNum=1; fileNum<=inFileCount; fileNum++)
			{
				// Determine the sample of combinations
				int allCountNum = choose(inFileCount, fileNum);
				int countNum = (allCountNum <= maxCom) ?allCountNum :maxCom;
				BitSet[] countSigs = new BitSet[countNum];
				for (int i=0; i<countNum; i++)
					countSigs[i] = new BitSet(inFileCount);

				if (allCountNum <= maxCom)	// Can count all of them
				{
					int[] countSig = new int[fileNum];
					for (int i=0; i<countNum; i++)
					{
						if (i == 0)
							for (int j=0; j<fileNum; j++)
								countSig[j] = j;
						else
						{
							// Determine the point of increment
							int incPt = fileNum - 1;
							while (countSig[incPt] == inFileCount + incPt - fileNum)
								incPt --;
							countSig[incPt] ++;
							for (int j=incPt+1; j<fileNum; j++)
								countSig[j] = countSig[j-1] + 1;
						}
						for (int j=0; j<fileNum; j++)
							countSigs[i].set(countSig[j]);
					}
				}
				else				// Can only count a sample of them
				{
					for (int i=0; i<countNum; i++)
					{
						int chosenCount = 0;
						while (chosenCount < fileNum)
						{
							int nextChosen = rand.nextInt(inFileCount);
							if (countSigs[i].get(nextChosen) == false)
							{
								countSigs[i].set(nextChosen);
								chosenCount ++;
							}
						}
					}
				}

				// Do the actual counting
				int[] ntCounts = new int[countNum];
				for (String id : fastEndPoints.keySet())
				{
					int[] fastIdEndPoints = fastEndPoints.get(id);
					BitSet[] coveredFeatures = coverMap.get(id);

					for (int i=0; i<fastIdEndPoints.length-1; i++)
					{
						int start = fastIdEndPoints[i];
						int endExc = fastIdEndPoints[i+1];
						int lengthInc = endExc - start;

						for (int j=0; j<countNum; j++)
							if (countSigs[j].intersects(coveredFeatures[i]))
								ntCounts[j] += lengthInc;
					}
				}

				// Store counts in lists
				java.util.List<Double> ntDist = new ArrayList<Double>(countNum);
				ntDists.add(ntDist);

				if (totNtCount == 0)
					for (int i=0; i<countNum; i++)
						ntDist.add(new Double(ntCounts[i]));
				else
					for (int i=0; i<countNum; i++)
						ntDist.add(new Double(ntCounts[i]/ totNtCount));
			}
		}

		// Output counts
		System.out.println("Outputing the count/fraction statistics.");
		PrintWriter pw = new PrintWriter(new FileWriter(outCountFile));
		pw.println("Number of features\t25 percentile\tMedian\t75 percentile");
		for (int i=0; i<inFileCount; i++)
		{
			java.util.List<Double> ntDist = ntDists.get(i);
			Collections.sort(ntDist);
			int countNum = ntDist.size();
			pw.println((i + 1) + "\t" +
			           ntDist.get((int)(countNum * 0.25)) + "\t" +
			           ntDist.get((int)(countNum * 0.50)) + "\t" +
			           ntDist.get((int)(countNum * 0.75)));
		}
		pw.close();

		// Convert data into JFreeChart format
		DefaultBoxAndWhiskerCategoryDataset ntData = (totNtCount == 0)
		                                           ? new DefaultBoxAndWhiskerCategoryDataset()
		                                           : new DefaultBoxAndWhiskerCategoryDatasetWithNoFarOuts();
		for (int i=0; i<inFileCount; i++)
		{
			String colID = "" + (i + 1);
			ntData.add(ntDists.get(i), "", colID);
		}

		// Create the chart
		System.out.println("Creating the plot.");
		JFreeChart chart = ChartFactory.createBoxAndWhiskerChart(
	        	null,
	        	"Number of features",
	        	(totNtCount == 0) ?"Number of covered nucleotides" :"Fraction of covered nucleotides",
	        	ntData,
	        	false);
	        CategoryPlot plot = (CategoryPlot)chart.getPlot();

		// Configure background color
		plot.setBackgroundPaint(Color.white);

		// Configure gridlines
		plot.setRangeGridlinesVisible(false);
	        plot.setDomainGridlinesVisible(false);

		// Set axis bounds
		if (totNtCount != 0)
		{
			plot.getRangeAxis().setLowerBound(0);
			plot.getRangeAxis().setUpperBound(1);
		}

		// Output document settings
		int leftMargin = 10;
		int rightMargin = 10;
		int topMargin = 10;
		int bottomMargin = 10;
		int documentWidth = 800;
		int documentHeight = 600;

		// Export chart to PDF
		FileOutputStream fos = new FileOutputStream(outGraphFile);
		com.lowagie.text.Rectangle rectangle = new com.lowagie.text.Rectangle(documentWidth, documentHeight);
		Document document = new Document(rectangle, leftMargin, rightMargin, topMargin, bottomMargin);
		PdfWriter pdfwriter = PdfWriter.getInstance(document, fos);
		document.open();
		PdfContentByte pdfContent = pdfwriter.getDirectContent();
		PdfTemplate pdfTemplate = pdfContent.createTemplate(documentWidth, documentHeight);
		Graphics2D g2d = pdfTemplate.createGraphics(documentWidth, documentHeight, new DefaultFontMapper());

		// Add the chart to the file
		Rectangle2D.Double offset = new Rectangle2D.Double(leftMargin, topMargin, documentWidth - leftMargin - rightMargin, documentHeight - topMargin - bottomMargin);
		chart.draw(g2d, offset);

		// Complete the plot
		g2d.dispose();
		pdfContent.addTemplate(pdfTemplate, 0, 0);
		document.close();
	}


	/*----------------------------------------------------------------------
	 Non-public methods
	----------------------------------------------------------------------*/

	/**
	 *	Get the number of ways to choose r balls from n balls.
	 *	@param		n			The number of balls in total
	 *	@param		r			The number of balls being chosen
	 *	@return					The number of ways
	 *	@exception	IllegalArgumentException	If n or r is negative,
	 *							or if n is smaller than r
	 */
	public static int choose(int n, int r) throws IllegalArgumentException
	{
		if (n < 0)
			throw new IllegalArgumentException("Negative n.");
		if (r < 0)
			throw new IllegalArgumentException("Negative r.");
		if (n < r)
			throw new IllegalArgumentException("n smaller than r.");

		double result = 1;
		int bound = (r > n - r) ?(n - r) :r;
		for (int i=0; i<bound; i++)
			result *= ((double)n - i) / (i + 1);
		return (int)result;
	}


	/*----------------------------------------------------------------------
	 Non-public classes
	----------------------------------------------------------------------*/

	protected static class DefaultBoxAndWhiskerCategoryDatasetWithNoFarOuts extends DefaultBoxAndWhiskerCategoryDataset
	{
		public Number getMinOutlier(int row, int column)
		{
			return 0;
		}

		public Number getMinOutlier(Comparable rowKey, Comparable columnKey)
		{
			return 0;
		}

		public Number getMaxOutlier(int row, int column)
		{
			return 1;
		}

		public Number getMaxOutlier(Comparable rowKey, Comparable columnKey)
		{
			return 1;
		}
    	}
}
