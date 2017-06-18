import java.io.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.print.*;
import javax.swing.*;

/**
 * ManageOutput
 <P>
 DataAnalysisJava is a first version of a data analysis package for 
 Java.  Class ManageOutput supports save and print functionality 
 from output display windows.
 <P>
 @author Fionn Murtagh, fmurtagh@acm.org
 @version 0.0, 2003-March-15
*/

public class ManageOutput extends JMenuBar {

public ManageOutput(final JTextArea outext) {

     JMenu fileMenu = new JMenu("File");

     ActionListener printListener = new ActionListener() {
	 public void  actionPerformed(ActionEvent event) {
                String cmd = event.getActionCommand();
                if (cmd.equals("Save...")) doSave(outext);
                if (cmd.equals("Print...")) doPrint(outext);
 	        System.out.println("Save/print completed...");
	 }
     };

     JMenuItem saveCmd = new JMenuItem("Save...");
     saveCmd.addActionListener(printListener);
     saveCmd.setAccelerator( KeyStroke.getKeyStroke("ctrl S") );
     fileMenu.add(saveCmd);

     JMenuItem printCmd = new JMenuItem("Print...");
     printCmd.addActionListener(printListener);
     printCmd.setAccelerator( KeyStroke.getKeyStroke("ctrl P") );
     fileMenu.add(printCmd);

     add(fileMenu);

}  // End of: ManageOutput contstructor

   private void doSave(JTextArea outext) {
          // Carry out the Save command by letting the user specify
          // an output file and writing the text from the JTextArea
          // to that file.
      File file;  // The file that the user wants to save.
      JFileChooser fd; // File dialog that lets the user specify the file.
      fd = new JFileChooser(new File("."));
      fd.setDialogTitle("Save Text As...");
      int action = fd.showSaveDialog(this);
      if (action != JFileChooser.APPROVE_OPTION) {
             // User has canceled, or an error occurred.
         return;
      }
      file = fd.getSelectedFile();
      if (file.exists()) {
            // If file already exists, ask before replacing it.
         action = JOptionPane.showConfirmDialog(this,
                     "Replace existing file?");
         if (action != JOptionPane.YES_OPTION)
            return;
      }
      try {
            // Create a PrintWriter for writing to the specified
            // file and write the text from the window to that stream.
         PrintWriter out = new PrintWriter(new FileWriter(file));
	 //JTextArea outext = DataAnalysis.outext1;
         String contents = outext.getText();
         out.print(contents);
         if (out.checkError())
            throw new IOException("Error while writing to file.");
         out.close();
      }
      catch (IOException e) {
            // Some error has occured while trying to write.
            // Show an error message.
         JOptionPane.showMessageDialog(this,
             "Sorry, an error has occurred:\n" + e.getMessage());
      }
   }  // End of: doSave



   private void doPrint(JTextArea outext) {
          // Carry out the Print command by letting the user specify
          // the printer and printing the text from the JTextArea.

       PrintUtilities.printComponent(outext);

    }  // End of: doPrint


}  // End of: ManageOutput class
