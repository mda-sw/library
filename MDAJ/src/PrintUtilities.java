import java.io.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.print.*;
import javax.swing.*;

/**
 * PrintUtilities
 <P>
 DataAnalysisJava is a first version of a data analysis package for 
 Java.  Class PrintUtilities supports print functionality 
 from output display windows.
 <P>
 @author Fionn Murtagh, fmurtagh@acm.org
 @version 0.0, 2003-March-15
*/

public class PrintUtilities implements Printable {
     private Component componentToBePrinted;

     public static void printComponent(Component c) {
         new PrintUtilities(c).print();
     }

     public PrintUtilities(Component componentToBePrinted) {
         this.componentToBePrinted = componentToBePrinted;
     }

     public void print() {
         PrinterJob printJob = PrinterJob.getPrinterJob();
	 PageFormat pf1 = printJob.defaultPage();
         printJob.setPrintable(this,pf1);
         if (printJob.printDialog()) try {
            printJob.print();
         } 
         catch(PrinterException pe) 
         {
            System.out.println("Error printing: " + pe);
         }
     }

     public int print(Graphics g, PageFormat pageFormat, int pageIndex) {
         if (pageIndex > 0) {
            return(NO_SUCH_PAGE);
         } else {
            Graphics2D g2d = (Graphics2D)g;
            g2d.translate(pageFormat.getImageableX(), 
                       pageFormat.getImageableY());
            disableDoubleBuffering(componentToBePrinted);
            componentToBePrinted.paint(g2d);
            enableDoubleBuffering(componentToBePrinted);
            return(PAGE_EXISTS);
	 }
     }

     /** The speed and quality of printing suffers dramatically if
      * any of the containers have double buffering turned on.
      * So this turns if off globally.
      * @see enableDoubleBuffering
      */
      public static void disableDoubleBuffering(Component c) {
          RepaintManager currentManager = RepaintManager.currentManager(c);
          currentManager.setDoubleBufferingEnabled(false);
      }

      /** Re-enables double buffering globally. */

     public static void enableDoubleBuffering(Component c) {
         RepaintManager currentManager = RepaintManager.currentManager(c);
         currentManager.setDoubleBufferingEnabled(true);
     }

}  // End of: PrintUtilities

