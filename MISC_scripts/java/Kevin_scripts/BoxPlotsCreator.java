import com.lowagie.text.*;
import com.lowagie.text.pdf.*;

import org.jfree.chart.*;
import org.jfree.chart.axis.*;
import org.jfree.chart.plot.*;
import org.jfree.chart.renderer.category.*;
import org.jfree.chart.title.*;
import org.jfree.data.statistics.*;
import org.jfree.ui.*;

import weka.core.*;
import weka.core.converters.*;

import java.awt.*;
import java.awt.geom.*;
import java.io.*;
import java.util.*;

/**
 *	<pre>
 *	This is a class that creates box plots from a Weka .arff file.
 *	</pre>
 *
 *	@author		Kevin Yuk-Lap Yip
 *	@version	1.0 (February 11, 2010)
 *
 *	<pre>
 *	Change History:
 *	1.0	- Initial version
 *	</pre>
 */

public class BoxPlotsCreator
{
	/*----------------------------------------------------------------------
	 Public methods
	----------------------------------------------------------------------*/

	public static void main(String[] argv) throws Exception
	{
		int argc = 0;
		File inFile	= new File(argv[argc++]);
		File outFile	= new File(argv[argc++]);

		// Read data
		ArffLoader dataLoader = new ArffLoader();
		dataLoader.setSource(inFile);
		Instances wData = dataLoader.getDataSet();
		wData.setClassIndex(wData.numAttributes() - 1);

		// Get the minimum and maximum values of each chosen attribute
		double[] minValues = new double[wData.numAttributes() - 6];
		double[] maxValues = new double[wData.numAttributes() - 6];

		for (int i=0; i<wData.numInstances(); i++)
		{
			Instance instance = wData.instance(i);
			for (int a=5; a<wData.numAttributes() - 1; a++)
			{
				double value = instance.value(a);
				if (value < minValues[a-5])
					minValues[a-5] = value;
				if (value > maxValues[a-5])
					maxValues[a-5] = value;
			}
		}

		double[] ranges = new double[wData.numAttributes() - 6];
		for (int a=0; a<wData.numAttributes() - 6; a++)
			ranges[a] = maxValues[a] - minValues[a];

		// Define class names
		String[] classNames = new String[wData.numClasses()];
		for (int c=0; c<wData.numClasses(); c++)
		{
			String className = wData.classAttribute().value(c);
			     if (className.equals("ncRNA_sampled"))		classNames[c] = "Known ncRNAs";
			else if (className.equals("exon_CDS"))			classNames[c] = "CDSs";
			else if (className.equals("UTR"))			classNames[c] = "UTRs";
			else if (className.equals("intergenic_gold"))		classNames[c] = "Intergenic Regions";

			else if (className.equals("miRNA"))			classNames[c] = "miRNAs";
			else if (className.equals("rRNA"))			classNames[c] = "rRNAs";
			else if (className.equals("scRNA"))			classNames[c] = "scRNAs";
			else if (className.equals("snoRNA"))			classNames[c] = "snoRNAs";
			else if (className.equals("snRNA"))			classNames[c] = "snRNAs";
			else if (className.equals("tRNA"))			classNames[c] = "tRNAs";

			else if (className.equals("Known_ncRNA"))		classNames[c] = "Known ncRNAs";
			else if (className.equals("Type_I_novel_ncRNA"))	classNames[c] = "High-confidence novel ncRNAs";
			else if (className.equals("Type_II_novel_ncRNA"))	classNames[c] = "Medium-confidence novel ncRNAs";

			else if (className.equals("ncRNA"))			classNames[c] = "Known ncRNAs";
			else if (className.equals("ncRNA_selected"))		classNames[c] = "Known ncRNAs";
			else if (className.equals("exon_CCDS"))			classNames[c] = "CCDSs";
			else if (className.equals("UTR"))			classNames[c] = "UTRs";
			else if (className.equals("ancestral_repeat"))		classNames[c] = "Ancestral repeats";

			else							classNames[c] = className;
		}

		// Define feature names
		String[] attrNames = new String[wData.numAttributes() - 6];
		for (int a=5; a<wData.numAttributes() - 1; a++)
		{
			String attrName = wData.attribute(a).name();
			     if (attrName.equals("gc"))				attrNames[a-5] = "GC percentage";
			else if (attrName.equals("identities"))			attrNames[a-5] = "DNA conservation";
			else if (attrName.equals("zscore_fix"))			attrNames[a-5] = "Predicted secondary structure free energy";
			else if (attrName.equals("sci_fix"))			attrNames[a-5] = "Predicted secondary structure conservation";
			else if (attrName.equals("tblastx_fix"))		attrNames[a-5] = "Predicted protein sequence conservation";
			else if (attrName.equals("polyA_RNAseq_max_all"))	attrNames[a-5] = "Poly-A+ RNAseq max.";
			else if (attrName.equals("small_RNAseq_max_all"))	attrNames[a-5] = "Small RNAseq max.";
			else if (attrName.equals("array_max_totalRNA"))		attrNames[a-5] = "Total RNA tiling array max.";
			else if (attrName.equals("array_max_polyA"))		attrNames[a-5] = "Poly-A+ RNA tiling array max.";

			else							attrNames[a-5] = attrName;
		}

		// Convert data into JFreeChart format
		DefaultBoxAndWhiskerCategoryDataset jData = new DefaultBoxAndWhiskerCategoryDatasetWithNoFarOuts();
		java.util.List[][] distributions = new java.util.List[wData.numClasses()][wData.numAttributes() - 6];
		for (int c=0; c<wData.numClasses(); c++)
			for (int a=5; a<wData.numAttributes() - 1; a++)
				distributions[c][a-5] = new ArrayList<Double>();

		for (int i=0; i<wData.numInstances(); i++)
		{
			Instance instance = wData.instance(i);
			if (instance.classIsMissing())
				continue;

			for (int a=5; a<wData.numAttributes() - 1; a++)
				distributions[(int)instance.classValue()][a-5].add(new Double((instance.value(a) - minValues[a-5]) / ranges[a-5]));
		}

		for (int c=0; c<wData.numClasses(); c++)
			for (int a=5; a<wData.numAttributes() - 1; a++)
				jData.add(distributions[c][a-5], classNames[c], attrNames[a-5]);

		// Create the chart
	        JFreeChart chart = ChartFactory.createBoxAndWhiskerChart("", "", "Normalized value", jData, true);
	        CategoryPlot plot = (CategoryPlot)chart.getPlot();

		// Configure background color
		plot.setBackgroundPaint(Color.white);

		// Configure legend location
		LegendTitle legends = chart.getLegend();
		legends.setPosition(RectangleEdge.TOP);

		// Configure gridlines
		plot.setRangeGridlinesVisible(false);
//	        plot.setRangePannable(true);
	        plot.setDomainGridlinesVisible(false);
//	        plot.setDomainGridlinePaint(Color.gray);
//	        plot.setDomainGridlinePosition(CategoryAnchor.END);

		// Configure range axis
	        NumberAxis numberAxis = (NumberAxis)plot.getRangeAxis();
	        numberAxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits());

		// Configure domain axis
	        CategoryAxis domainAxis = (CategoryAxis)plot.getDomainAxis();
		domainAxis.setLowerMargin(0.02);
		domainAxis.setUpperMargin(0.02);
		domainAxis.setCategoryMargin(0.20);
		for (int a=1; a<attrNames.length; a+=2)
		{
			CategoryMarker marker = new CategoryMarker(attrNames[a]);
			marker.setPaint(new Color(203, 203, 203));
			plot.addDomainMarker(marker, Layer.BACKGROUND);
		}

		domainAxis.setMaximumCategoryLabelLines(5);
		domainAxis.setMaximumCategoryLabelWidthRatio(1.3f);
	        for (int a=0; a<attrNames.length; a++)
		        domainAxis.setTickLabelFont(attrNames[a], new java.awt.Font("Arial", java.awt.Font.PLAIN, 11));

		// Configure series colors
		BoxAndWhiskerRenderer renderer = (BoxAndWhiskerRenderer) plot.getRenderer();
		for (int c=0; c<classNames.length; c++)
		{
			Color seriesColor = null;

			     if (classNames[c].equals("Known ncRNAs"))			seriesColor = new Color( 63,  63, 191);
			else if (classNames[c].equals("CDSs"))				seriesColor = new Color(223,  95,  95);
			else if (classNames[c].equals("UTRs"))				seriesColor = new Color(223, 223,  95);
			else if (classNames[c].equals("Intergenic Regions"))		seriesColor = new Color( 95, 223,  95);

			else if (classNames[c].equals("miRNAs"))			seriesColor = new Color(160,  33, 240);
			else if (classNames[c].equals("rRNAs"))				seriesColor = new Color(223,  95,  95);
			else if (classNames[c].equals("scRNAs"))			seriesColor = new Color( 95, 223,  95);
			else if (classNames[c].equals("snoRNAs"))			seriesColor = new Color( 95, 223, 223);
			else if (classNames[c].equals("snRNAs"))			seriesColor = new Color(223, 223,  95);
			else if (classNames[c].equals("tRNAs"))				seriesColor = new Color( 95,  95, 223);

			else if (classNames[c].equals("Known ncRNAs"))			seriesColor = new Color( 63,  63, 191);
			else if (classNames[c].equals("High-confidence novel ncRNAs"))	seriesColor = new Color( 95, 223, 223);
			else if (classNames[c].equals("Medium-confidence novel ncRNAs"))seriesColor = new Color(223,  95, 223);

			else if (classNames[c].equals("Known ncRNAs"))			seriesColor = new Color( 63,  63, 191);
			else if (classNames[c].equals("CCDSs"))				seriesColor = new Color(223,  95,  95);
			else if (classNames[c].equals("UTRs"))				seriesColor = new Color(223, 223,  95);
			else if (classNames[c].equals("Ancestral repeats"))		seriesColor = new Color( 95, 223,  95);

			renderer.setSeriesPaint(c, seriesColor);
		}

		// Output document settings
		int leftMargin = 10;
		int rightMargin = 10;
		int topMargin = 10;
		int bottomMargin = 10;
//		int documentWidth = 800;
		int documentWidth = 2400;
		int documentHeight = 400;

		// Export chart to PDF
		FileOutputStream fos = new FileOutputStream(outFile);
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
