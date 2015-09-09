import com.lowagie.text.*;
import com.lowagie.text.pdf.*;

import org.jfree.chart.*;
import org.jfree.chart.plot.*;
import org.jfree.chart.renderer.xy.*;
import org.jfree.chart.title.*;
import org.jfree.data.xy.*;
import org.jfree.ui.*;

import weka.core.*;
import weka.core.converters.*;

import java.awt.*;
import java.awt.geom.*;
import java.io.*;
import java.util.*;

/**
 *	<pre>
 *	This is a class that creates scatter plots.
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

public class ScatterPlotsCreator
{
	/*----------------------------------------------------------------------
	 Public methods
	----------------------------------------------------------------------*/

	public static void main(String[] argv) throws Exception
	{
		int argc = 0;
		File inFile		= new File(argv[argc++]);
		String coorIndicesStr	= argv[argc++];
		double jitter		= Double.parseDouble(argv[argc++]);
		File outFile		= new File(argv[argc++]);

		StringTokenizer st = new StringTokenizer(coorIndicesStr, ",", false);
		int[] coorIndices = new int[st.countTokens()];
		for (int i=0; i<coorIndices.length; i++)
			coorIndices[i] = Integer.parseInt(st.nextToken());

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
			else if (className.equals("intergenic_gold"))
			{
				if (outFile.getName().contains("Fig4Aleft.pdf"))	classNames[c] = "Unepxressed Intergenic";
				else							classNames[c] = "Intergenic Regions";
			}

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

			else if (attrName.equals("ncRNA_score"))		attrNames[a-5] = "ncRNA score";
			else if (attrName.equals("CDS_score"))			attrNames[a-5] = "CDS score";

			else							attrNames[a-5] = attrName;
		}

		// Count the number of instances in each class
		int[] classInstanceCounts = new int[wData.numClasses()];
		for (int i=0; i<wData.numInstances(); i++)
		{
			Instance instance = wData.instance(i);
			if (!instance.classIsMissing())
				classInstanceCounts[(int)instance.classValue()] ++;
		}

		// Convert data into JFreeChart format
		DefaultXYDataset jData = new DefaultXYDataset();
		double[][][] distributions = new double[wData.numClasses()][coorIndices.length][];
		for (int c=0; c<wData.numClasses(); c++)
			for (int a=0; a<coorIndices.length; a++)
				distributions[c][a] = new double[classInstanceCounts[c]];

		int[] classPtrs = new int[wData.numClasses()];
		Random rand = new Random(0);
		for (int i=0; i<wData.numInstances(); i++)
		{
			Instance instance = wData.instance(i);
			if (instance.classIsMissing())
				continue;
			int classValue = (int)instance.classValue();

			for (int a=0; a<coorIndices.length; a++)
			{
				double value = instance.value(coorIndices[a]);
				if (jitter != 0)
				{
					value += (rand.nextDouble() - 0.5) * jitter;
					if (value < 0)
						value = 0;
					if (outFile.getName().contains("Fig4A"))
					{
						if (a == 0 && value > 0.18 && instance.value(coorIndices[a]) <= 0.18)
							value = instance.value(coorIndices[a]);
						if (a == 0 && value < 0.69 && instance.value(coorIndices[a]) >= 0.69)
							value = instance.value(coorIndices[a]);
					}
					else if (outFile.getName().contains("human.pdf"))
					{
						if (a == 0 && value > 0.07 && instance.value(coorIndices[a]) <= 0.07)
							value = instance.value(coorIndices[a]);
						if (a == 0 && value < 0.68 && instance.value(coorIndices[a]) >= 0.68)
							value = instance.value(coorIndices[a]);
					}
					else if (outFile.getName().contains("prediction_gold_gold_scatterplots.pdf"))
					{
						if (a == 0 && value > 0.21 && instance.value(coorIndices[a]) <= 0.21)
							value = instance.value(coorIndices[a]);
						if (a == 0 && value < 0.68 && instance.value(coorIndices[a]) >= 0.68)
							value = instance.value(coorIndices[a]);
					}
				}
				distributions[classValue][a][classPtrs[classValue]] = value;
			}
			classPtrs[classValue] ++;
		}

		for (int c=0; c<wData.numClasses(); c++)
			jData.addSeries(classNames[c], distributions[c]);

		// Create the chart
		boolean showLegend = (outFile.getName().contains("Fig4A") || outFile.getName().contains("human.pdf") || outFile.getName().contains("prediction_gold_gold_scatterplots.pdf"));
		JFreeChart chart = ChartFactory.createScatterPlot("", attrNames[coorIndices[0]-5], attrNames[coorIndices[1]-5], jData, PlotOrientation.VERTICAL, showLegend, false, false);
	        XYPlot plot = (XYPlot)chart.getPlot();

		// Configure background color
		plot.setBackgroundPaint(Color.white);

		// Configure gridlines
		plot.setRangeGridlinesVisible(false);
	        plot.setDomainGridlinesVisible(false);

		// Configure axis fonts
		java.awt.Font axesFont = new java.awt.Font("Arial", java.awt.Font.PLAIN, 8);
		plot.getDomainAxis().setTickLabelFont(axesFont);
		plot.getRangeAxis().setTickLabelFont(axesFont);

		// Configure series colors
		XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) plot.getRenderer();
		for (int c=0; c<classNames.length; c++)
		{
			Color seriesColor = null;

			     if (classNames[c].equals("Known ncRNAs"))			seriesColor = new Color( 63,  63, 191);
			else if (classNames[c].equals("CDSs"))				seriesColor = new Color(223,  95,  95);
			else if (classNames[c].equals("UTRs"))				seriesColor = new Color(223, 223,  95);
			else if (classNames[c].equals("Intergenic Regions"))		seriesColor = new Color( 95, 223,  95);
			else if (classNames[c].equals("Unepxressed Intergenic"))	seriesColor = new Color( 95, 223,  95);

			else if (classNames[c].equals("miRNAs"))			seriesColor = new Color(160,  33, 240);
			else if (classNames[c].equals("rRNAs"))				seriesColor = new Color(223,  95,  95);
			else if (classNames[c].equals("scRNAs"))			seriesColor = new Color( 95, 223,  95);
			else if (classNames[c].equals("snoRNAs"))			seriesColor = new Color( 95, 223, 223);
			else if (classNames[c].equals("snRNAs"))			seriesColor = new Color(223, 223,  95);
			else if (classNames[c].equals("tRNAs"))				seriesColor = new Color( 95,  95, 223);

			else if (classNames[c].equals("Known ncRNAs"))			seriesColor = new Color( 63,  63, 191);
			else if (classNames[c].equals("High-confidence novel ncRNAs"))
			{
				if (outFile.getName().contains("grayscale"))		seriesColor = new Color( 63,  63,  63);
				else							seriesColor = new Color( 95, 223, 223);
			}
			else if (classNames[c].equals("Medium-confidence novel ncRNAs"))
			{
				if (outFile.getName().contains("grayscale"))		seriesColor = new Color(191, 191, 191);
				else							seriesColor = new Color(223,  95, 223);
			}

			else if (classNames[c].equals("Known ncRNAs"))			seriesColor = new Color( 63,  63, 191);
			else if (classNames[c].equals("CCDSs"))				seriesColor = new Color(223,  95,  95);
			else if (classNames[c].equals("UTRs"))				seriesColor = new Color(223, 223,  95);
			else if (classNames[c].equals("Ancestral repeats"))		seriesColor = new Color( 95, 223,  95);

			renderer.setSeriesPaint(c, seriesColor);
		}

		// Configure series shapes, sizes and fills
		for (int c=0; c<classNames.length; c++)
		{
			renderer.setSeriesShape(c, new Ellipse2D.Double(0, 0, 1, 1));
			if (c == 0 && classNames[c].equals("Known ncRNAs") && !classNames[1].contains("Novel") &&
			    !outFile.getName().contains("Fig4A") &&
			    !outFile.getName().contains("human.pdf") &&
			    !outFile.getName().contains("prediction_gold_gold_scatterplots.pdf"))
				renderer.setSeriesShape(c, new Ellipse2D.Double(0, 0, 2, 2));
			renderer.setSeriesShapesFilled(c, true);
		}

		if (outFile.getName().contains("Fig4A") ||
		    outFile.getName().contains("human.pdf") ||
		    outFile.getName().contains("prediction_gold_gold_scatterplots.pdf"))
		{
			// Set axis bounds
			plot.getDomainAxis().setLowerBound(0);
			plot.getDomainAxis().setUpperBound(1);
			plot.getRangeAxis().setLowerBound(0);
			plot.getRangeAxis().setUpperBound(1);

			// Set legends
			LegendItemCollection oldLegendItems = plot.getLegendItems();
			LegendItemCollection newLegendItems = new LegendItemCollection();
			for (int i=0; i<oldLegendItems.getItemCount(); i++)
			{
				LegendItem oldLegendItem = oldLegendItems.get(i);
				if (outFile.getName().contains("Fig4Aleft.pdf") ||
				    outFile.getName().contains("human.pdf") ||
				    outFile.getName().contains("prediction_gold_gold_scatterplots.pdf"))
				{
					if (i == 2)	oldLegendItem = oldLegendItems.get(3);
					if (i == 3)	oldLegendItem = oldLegendItems.get(2);
				}

				LegendItem newLegendItem = new LegendItem(oldLegendItem.getLabel(),
				                                          oldLegendItem.getDescription(),
				                                          oldLegendItem.getToolTipText(),
				                                          oldLegendItem.getURLText(),
				                                          new Rectangle2D.Double(0, 0, 10, 10),
				                                          oldLegendItem.getFillPaint(),
				                                          oldLegendItem.getOutlineStroke(),
				                                          oldLegendItem.getOutlinePaint());
				newLegendItems.add(newLegendItem);
			}
			LegendTitle legends = chart.getLegend();
			legends.setPosition(RectangleEdge.BOTTOM);
			legends.setItemFont(new java.awt.Font("Arial", java.awt.Font.PLAIN, 7));
			plot.setFixedLegendItems(newLegendItems);

			// Add markers
			Color markerColor = new Color(63, 63, 63);
			Stroke markerStroke = new BasicStroke(1, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 10.0f, new float[] {10.0f, 10.0f}, 0.0f);

			if (outFile.getName().contains("Fig4A"))
			{
				ValueMarker lowMarker = new ValueMarker(0.18);
				lowMarker.setPaint(markerColor);
				lowMarker.setStroke(markerStroke);
				plot.addDomainMarker(lowMarker, Layer.BACKGROUND);

				ValueMarker highMarker = new ValueMarker(0.69);
				highMarker.setPaint(markerColor);
				highMarker.setStroke(markerStroke);
				plot.addDomainMarker(highMarker, Layer.BACKGROUND);
			}
			else if (outFile.getName().contains("human.pdf"))
			{
				ValueMarker lowMarker = new ValueMarker(0.07);
				lowMarker.setPaint(markerColor);
				lowMarker.setStroke(markerStroke);
				plot.addDomainMarker(lowMarker, Layer.BACKGROUND);

				ValueMarker highMarker = new ValueMarker(0.68);
				highMarker.setPaint(markerColor);
				highMarker.setStroke(markerStroke);
				plot.addDomainMarker(highMarker, Layer.BACKGROUND);
			}
			else if (outFile.getName().contains("prediction_gold_gold_scatterplots.pdf"))
			{
				ValueMarker lowMarker = new ValueMarker(0.21);
				lowMarker.setPaint(markerColor);
				lowMarker.setStroke(markerStroke);
				plot.addDomainMarker(lowMarker, Layer.BACKGROUND);

				ValueMarker highMarker = new ValueMarker(0.68);
				highMarker.setPaint(markerColor);
				highMarker.setStroke(markerStroke);
				plot.addDomainMarker(highMarker, Layer.BACKGROUND);
			}
		}

		// Output document settings
		int leftMargin = 10;
		int rightMargin = 10;
		int topMargin = 10;
		int bottomMargin = 10;
		int documentWidth = 400;
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
}
