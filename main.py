# Imports
from Bio.Align.Applications import ClustalwCommandline  # For MSA
from mainWindow import Ui_MainWindow  # GUI file
from PyQt5 import QtWidgets  # For GUI
from Bio.SeqUtils import GC  # To calculate GC Content
from Bio import Phylo  # For Drawing Phylognetic Tree
from Bio import SeqIO  # For reading Fasta Files
from Bio.Seq import *
import sys
import os
from drawFunc import draw  # For Drawing Phylogentic Tree on the GUI
import cv2  # For Reading an Image
import pyqtgraph as pg  # For Reading an Image as a graph on the GUI
from Bio import pairwise2  # For Pairwise Sequence Alignment
from Bio.pairwise2 import format_alignment  # For formating alignment


class ApplicationWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super(ApplicationWindow, self).__init__()
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.filepath1 = ""
        self.filepath2 = ""
        self.filepath3 = ""
        self.filepath4 = ""
        self.filepath5 = ""
        # Tab 1 buttons
        self.ui.browseBtn.clicked.connect(self.browseFile1)
        self.ui.resultBtn.clicked.connect(
            lambda: self.showResult1(filepath=self.filepath1))
        # Tab 2 buttons
        self.ui.extractBtn.hide()
        self.ui.browseBtn_2.clicked.connect(self.browseFile2)
        self.ui.resultBtn_2.clicked.connect(
            lambda: self.showResult2(filepath=self.filepath2))
        self.ui.extractBtn.clicked.connect(
            lambda: self.showResult2(filepath=self.filepath2, extract=1))
        # Tab 3 buttons
        self.ui.browseBtn_3.clicked.connect(self.browseFile3)
        self.ui.resultBtn_3.clicked.connect(
            lambda: self.showResult3(filepath=self.filepath3))
        # Tab 4 buttons
        self.ui.browseBtn_4.clicked.connect(self.browseFile4)
        self.ui.browseBtn_5.clicked.connect(self.browseFile5)
        self.ui.resultBtn_4.clicked.connect(lambda: self.showResult4(
            filepath1=self.filepath4, filepath2=self.filepath5))

    # ------------------------------------------------------------------------------------------------------------------ #
    # ----------------------------- TAB 3 (Calculate % of Chemical Components in DNA or RNA Sequence ------------------- #
    # ------------------------------------------------------------------------------------------------------------------ #
    def browseFile1(self):  # Browse function: browse file and check if file is fasta file or not
        self.filepath1 = QtWidgets.QFileDialog.getOpenFileName(
            self, 'Open file')
        if self.filepath1[0] != "":
            filename, extension = os.path.splitext(self.filepath1[0])
            if extension == ".fasta" or extension == ".FASTA" or extension == ".Fasta":
                self.ui.label.setText(self.filepath1[0])
            else:
                choice = QtWidgets.QMessageBox.question(
                    self, 'WARNING!', "Please Choose Fasta File", QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
                if choice == QtWidgets.QMessageBox.Yes:
                    self.browseTab1()
                elif choice == QtWidgets.QMessageBox.No:
                    return

    # Show Result: Calculations for Chemical Components of DNA & RNA
    def showResult1(self, filepath):
        if filepath != "":
            filepath = filepath[0]
        self.ui.outputTextView.clear()
        # Get current text in ui dropdown menu
        checkSequenceType = self.ui.comboBox.currentText()
        inputVal = self.ui.inputTextView.toPlainText()
        # Check if user choose one of Sequences type (DNA, RNA, Protein)
        if checkSequenceType == "DNA" or checkSequenceType == "RNA" or checkSequenceType == "Protein":
            if filepath != "" and inputVal == "":
                record = SeqIO.parse(filepath, 'fasta')
                for element in record:
                    aConent = element.seq.count('A')  # Count of (A) Content
                    cContent = element.seq.count('C')  # Count of (C) Content
                    gContent = element.seq.count('G')  # Count of (G) Content
                    seqLength = len(element.seq)
                    self.ui.outputTextView.append(
                        "Sequence: " + element.name)
                    self.ui.outputTextView.append(
                        "Description: " + element.description)
                    self.ui.outputTextView.append(
                        "Sequence Length: " + str(seqLength))
                    self.ui.outputTextView.append(
                        "Percentage of (A): " + str((aConent / seqLength) * 100))  # Calculates % of (A) Content
                    self.ui.outputTextView.append(
                        "Percentage of (C): " + str((cContent / seqLength) * 100))  # Calculates % of (C) Content
                    self.ui.outputTextView.append(
                        "Percentage of (G): " + str((gContent / seqLength) * 100))  # Calculates % of (G) Content
                    # Check if DNA is Chosen, then:
                    if checkSequenceType == "DNA":
                        tContent = element.seq.count(
                            'T')  # Count of (T) Content
                        self.ui.outputTextView.append(
                            "Percentage of (T): " + str((tContent / seqLength) * 100))  # Calculates % of (T) Content
                    # Check if RNA is Chosen, then:
                    if checkSequenceType == "RNA":
                        uContent = element.seq.count(
                            'U')  # Count of (U) Content
                        self.ui.outputTextView.append(
                            "Percentage of (U): " + str((uContent / seqLength) * 100))
                    self.ui.outputTextView.append(
                        "Percentage of (CG) Content: " + str(GC(element.seq)))  # Calculates % of (GC) Content
                    self.ui.outputTextView.append(
                        "--------------------------------------------------")
            # Check if input is text
            if inputVal != "" or (inputVal != "" and filepath != ""):
                self.ui.label.clear()
                if checkSequenceType == "DNA":  # Check if type of Sequence is DNA
                    DNASeq = "ATGC"
                    for i in inputVal:
                        if i in DNASeq.lower() or i in DNASeq.upper():  # Check in input text contains (A, C, G, T) only
                            self.ui.outputTextView.clear()
                            seq = inputVal
                            seq = seq.upper()
                            aConent = seq.count('A')  # Count of (A) Content
                            cContent = seq.count('C')  # Count of (C) Content
                            gContent = seq.count('G')  # Count of (G) Content
                            tContent = seq.count('T')  # Count of (T) Content
                            seqLength = len(seq)
                            self.ui.outputTextView.append(
                                "Description: " + "DNA Sequence")
                            self.ui.outputTextView.append(
                                "Sequence Length: " + str(seqLength))
                            self.ui.outputTextView.append(
                                "Percentage of (A): " + str((aConent / seqLength) * 100))  # Calculates % of (A) Content
                            self.ui.outputTextView.append(
                                "Percentage of (C): " + str((cContent / seqLength) * 100))  # Calculates % of (C) Content
                            self.ui.outputTextView.append(
                                "Percentage of (G): " + str((gContent / seqLength) * 100))  # Calculates % of (G) Content
                            self.ui.outputTextView.append(
                                "Percentage of (T): " + str((tContent / seqLength) * 100))  # Calculates % of (T) Content
                            self.ui.outputTextView.append(
                                "Percentage of (CG) Content: " + str(GC(seq)))  # Calculates % of (GC) Content
                            self.ui.outputTextView.append(
                                "--------------------------------------------------")
                        else:
                            self.ui.outputTextView.clear()
                            self.ui.outputTextView.setText("ERROR!")
                            break

                if checkSequenceType == "RNA":  # Check if type of Sequence is RNA
                    # Check in input text contains (A, C, G, U) only
                    RNASeq = "AUGC"
                    for i in inputVal:
                        if i in RNASeq.lower() or i in RNASeq.upper():
                            self.ui.outputTextView.clear()
                            seq = inputVal
                            seq = seq.upper()
                            aConent = seq.count('A')  # Count of (A) Content
                            cContent = seq.count('C')  # Count of (C) Content
                            gContent = seq.count('G')  # Count of (G) Content
                            uContent = seq.count('U')  # Count of (U) Content
                            seqLength = len(seq)
                            self.ui.outputTextView.append(
                                "Description: " + "RNA Sequence")
                            self.ui.outputTextView.append(
                                "Sequence Length: " + str(seqLength))
                            self.ui.outputTextView.append(
                                "Percentage of (A): " + str((aConent / seqLength) * 100))  # Calculates % of (A) Content
                            self.ui.outputTextView.append(
                                "Percentage of (C): " + str((cContent / seqLength) * 100))  # Calculates % of (C) Content
                            self.ui.outputTextView.append(
                                "Percentage of (G): " + str((gContent / seqLength) * 100))  # Calculates % of (G) Content
                            self.ui.outputTextView.append(
                                "Percentage of (U): " + str((uContent / seqLength) * 100))  # Calculates % of (U) Content
                            self.ui.outputTextView.append(
                                "Percentage of (CG) Content: " + str(GC(seq)))  # Calculates % of (CG) Content
                            self.ui.outputTextView.append(
                                "--------------------------------------------------")
                        else:
                            self.ui.outputTextView.clear()
                            self.ui.outputTextView.setText("ERROR!")
                            break
        else:
            choice = QtWidgets.QMessageBox.question(
                self, 'WARNING!', "Please Choose Method first", QtWidgets.QMessageBox.Ok)
            if choice == QtWidgets.QMessageBox.Ok:
                return
    # ------------------------------------------------------------------------------------------------------------------ #
    # -------------------------------------------------- End of TAB 3 -------------------------------------------------- #
    # ------------------------------------------------------------------------------------------------------------------ #

    # ------------------------------------------------------------------------------------------------------------------- #
    # ------------------------------------------ TAB 4 (Transcribe, Translate) ------------------------------------------ #
    # ------------------------------------------------------------------------------------------------------------------- #

    def browseFile2(self):  # Browse function: browse file and check if file is fasta file or not
        self.filepath2 = QtWidgets.QFileDialog.getOpenFileName(
            self, 'Open file')
        if self.filepath2[0] != "":
            filename, extension = os.path.splitext(self.filepath2[0])
            if extension == ".fasta" or extension == ".FASTA" or extension == ".Fasta":
                self.ui.label_2.setText(self.filepath2[0])
            else:
                choice = QtWidgets.QMessageBox.question(
                    self, 'WARNING!', "Please Choose Fasta File", QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
                if choice == QtWidgets.QMessageBox.Yes:
                    self.browseFile2()
                elif choice == QtWidgets.QMessageBox.No:
                    return

    # Show Result: Conversion of DNA Sequence to mRNA and Protein
    def showResult2(self, filepath, extract=-1):
        if filepath != "":
            filepath = filepath[0]
        self.ui.outputTextView_2.clear()
        checkMethod = self.ui.comboBox_2.currentText()
        inputVal = self.ui.inputTextView_2.toPlainText()
        # Check if user Chose a Transcribe / Translate methods from a drop down in the GUI
        if checkMethod == "Transcribe" or checkMethod == "Translate":
            if filepath != "" and inputVal == "":  # If input is a FASTA file
                record = SeqIO.parse(filepath, 'fasta')
                for element in record:
                    # If user chose "Transcribe", then:
                    if checkMethod == "Transcribe":
                        # Transcribe DNA Sequence
                        rnaSeq = transcribe(str(element.seq))
                        # Output of will appear on the GUI Output Text
                        self.ui.outputTextView_2.append(
                            "Sequence: " + element.name)  # Display Sequence Name/ID
                        self.ui.outputTextView_2.append(
                            "Description: " + element.description)  # Display Sequence Description
                        self.ui.outputTextView_2.append(
                            "RNA Sequence: " + rnaSeq)  # Display Sequence mRNA (Transcribed Sequence)
                        self.ui.outputTextView_2.append(
                            "--------------------------------------------------")
                        if extract == 1:
                            self.extractSequences(
                                filename=filepath, outputFile="rna_output.fasta")  # Generate Output Fasta file that conatins (Transcribed DNA)

                    # If user chose "Translate", then:
                    elif checkMethod == "Translate":
                        # Translate DNA Sequence
                        proteinSeq = translate(str(element.seq))
                        # Output of will appear on the GUI Output Text
                        self.ui.outputTextView_2.append(
                            "Sequence: " + element.name)  # Display Sequence Name/ID
                        self.ui.outputTextView_2.append(
                            "Description: " + element.description)  # Display Sequence Description
                        self.ui.outputTextView_2.append(
                            "RNA Sequence: " + proteinSeq)  # Display Sequence Protein (Translated Sequence)
                        self.ui.outputTextView_2.append(
                            "--------------------------------------------------")
                        if extract == 1:
                            self.extractSequences(
                                filename=filepath, outputFile="protein_output.fasta")  # Generate Output Fasta file that conatins (Transcribed DNA)

                self.ui.extractBtn.setVisible(True)

             # If input is a text
            if inputVal != "" or (inputVal != "" and filepath != ""):
                self.ui.label_2.clear()
                seq = inputVal
                seq = seq.upper()
                # If user chose "Transcribe", then:
                if checkMethod == "Transcribe":
                    # Transcribe DNA Sequence
                    rnaSeq = transcribe(str(seq))
                    # Output of will appear on the GUI Output Text
                    self.ui.outputTextView_2.append(
                        "Transcribed Sequence: " + rnaSeq)  # Display Sequence mRNA (Transcribed Sequence)
                    self.ui.outputTextView_2.append(
                        "--------------------------------------------------")

                # If user chose "Translate", then:
                elif checkMethod == "Translate":
                    # Translate DNA Sequence
                    proteinSeq = translate(str(seq))
                    # Output of will appear on the GUI Output Text
                    self.ui.outputTextView_2.append(
                        "Translated Sequence: " + proteinSeq)  # Display Sequence Protein (Translated Sequence)
                    self.ui.outputTextView_2.append(
                        "--------------------------------------------------")

        # Warn the User if the user doesn't Choose a method (Transcribe / Translate)
        else:
            choice = QtWidgets.QMessageBox.question(
                self, 'WARNING!', "Please Choose Method first", QtWidgets.QMessageBox.Ok)
            if choice == QtWidgets.QMessageBox.Ok:
                return

        if extract == 1:
            QtWidgets.QMessageBox.information(
                self, "Done", "File Generated Successfully :)")  # Message appear when output file extracted

    # Function that Extract a Transcribe or Translate Result file in the same directory of run file
    def extractSequences(self,  filename, outputFile):
        checkSequenceType = self.ui.comboBox_2.currentText()
        file = open(filename)
        fout = open(outputFile, 'w')
        record = SeqIO.parse(filename, 'fasta')
        for element in record:
            if checkSequenceType == "Transcribe":
                rnaSeq = transcribe(str(element.seq))
                for index, record in enumerate(SeqIO.parse(file, 'fasta')):
                    fout.write(">" + record.id + "\n")
                    fout.write(rnaSeq + "\n")

            elif checkSequenceType == "Translate":
                proteinSeq = translate(str(element.seq))
                for index, record in enumerate(SeqIO.parse(file, 'fasta')):
                    fout.write(">" + record.id + "\n")
                    fout.write(proteinSeq + "\n")
    # ------------------------------------------------------------------------------------------------------------------- #
    # --------------------------------------------------- End ofTAB 4 --------------------------------------------------- #
    # ------------------------------------------------------------------------------------------------------------------- #

    # ------------------------------------------------------------------------------------------------------------------- #
    # -------------------------------------- TAB 2 => MSA (Alignment / Phylo Tree) -------------------------------------- #
    # ------------------------------------------------------------------------------------------------------------------- #
    def browseFile3(self): # Browse function: browse file and check if file is fasta file or not
        self.filepath3 = QtWidgets.QFileDialog.getOpenFileName(
            self, 'Open file')
        if self.filepath3[0] != "":
            filename, extension = os.path.splitext(self.filepath3[0])
            if extension == ".fasta" or extension == ".FASTA" or extension == ".Fasta":
                self.ui.filepathView_3.setText(self.filepath3[0])
            else:
                choice = QtWidgets.QMessageBox.question(
                    self, 'WARNING!', "Please Choose Fasta File", QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
                if choice == QtWidgets.QMessageBox.Yes:
                    self.browseFile3()
                elif choice == QtWidgets.QMessageBox.No:
                    return

    # Show Result: Multiple Sequence Alignmen & Phylogenetic Tree
    def showResult3(self, filepath):
        clustalw_exe = r"C:\Program Files (x86)\ClustalW2\clustalw2.exe" # Clustalw2 path
        if filepath != "":
            filepath = filepath[0]
        self.ui.msaOutputText.clear()
        # Get Current Method from a Dropdown menu (in GUI)
        checkMethod = self.ui.comboBox_3.currentText()
        filename, extension = os.path.splitext(filepath)
        MSA_outputfile = "MSA_outalign.fasta" # Output Multiple Sequence Alignment file result
        phylo_outputfile = filename + ".ph" # Output Phylogentic file result
        # Check if user Chose a Sequence Alignment / Phylogenetic Tree methods from a drop down in the GUI 
        if checkMethod == "Sequence Alignment" or checkMethod == "Phylogenetic Tree":
            # If user chose "Sequence Alignment", then:
            if checkMethod == "Sequence Alignment":
                self.ui.graphicsView.clear()
                # MSA using Clustal2 and Generating output FASTA file of result
                cmd = ClustalwCommandline(clustalw_exe, infile=filepath,
                                          type="DNA", output='FASTA', outfile=MSA_outputfile)
                std_out, std_err = cmd()

                record = SeqIO.parse(MSA_outputfile, 'fasta')
                for element in record:
                    # Output of will appear on the GUI Output Text
                    self.ui.msaOutputText.append(
                        "Sequence ID: " + str(element.name)) # Display a Name of Sequence
                    self.ui.msaOutputText.append(
                        "Sequence: " + str(element.seq)) # Display a result of MSA 
                    self.ui.msaOutputText.append(
                        "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

            if checkMethod == "Phylogenetic Tree":
                # Phylogentic Tree using Clustal2 and Generating output file of result
                cmd = ClustalwCommandline(clustalw_exe, infile=filepath,
                                          type="DNA", clustering="NJ", tree=True, outputtree='phylip', output='FASTA')
                std_out, std_err = cmd()
                # Reading Ouput Tree File
                readTree = Phylo.read(phylo_outputfile, 'newick')
                # Draw the Tree on The GUI 
                # Draw Function is imported from "drawFunc.py": It's exactly a function in the documentation of Biopython but we edited a function to draw a picture of result on th ui
                draw(readTree)
                imgArr = cv2.imread("tree.jpg")
                img = pg.ImageItem(imgArr)
                img.rotate(270)
                self.ui.graphicsView.addItem(img) # Displaying an image of Phylogentic Tree
                self.ui.msaOutputText.append(std_out) # Displaying information about Output tree on a output text are in the GUI
        else:
            choice = QtWidgets.QMessageBox.question(
                self, 'WARNING!', "Please Choose Method first", QtWidgets.QMessageBox.Ok)
            if choice == QtWidgets.QMessageBox.Ok:
                return
    # ------------------------------------------------------------------------------------------------------------------ #
    # -------------------------------------------------- End of TAB 2 -------------------------------------------------- #
    # ------------------------------------------------------------------------------------------------------------------ #



    # ------------------------------------------------------------------------------------------------------------------ #
    # ---------------------------------- TAB 1 => Pairwise Alignemnt (Global - Local) ---------------------------------- #
    # ------------------------------------------------------------------------------------------------------------------ #
    def browseFile4(self): # Browse function: browse input file 1 and check if file is fasta file or not
        self.filepath4 = QtWidgets.QFileDialog.getOpenFileName(
            self, 'Open file')
        if self.filepath4[0] != "":
            filename, extension = os.path.splitext(self.filepath4[0])
            if extension == ".fasta" or extension == ".FASTA" or extension == ".Fasta":
                self.ui.filepathView_4.setText(self.filepath4[0])
            else:
                choice = QtWidgets.QMessageBox.question(
                    self, 'WARNING!', "Please Choose Fasta File", QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
                if choice == QtWidgets.QMessageBox.Yes:
                    self.browseFile4()
                elif choice == QtWidgets.QMessageBox.No:
                    return

    def browseFile5(self): # Browse function: browse input file 2 and check if file is fasta file or not
        self.filepath5 = QtWidgets.QFileDialog.getOpenFileName(
            self, 'Open file')
        if self.filepath5[0] != "":
            filename, extension = os.path.splitext(self.filepath5[0])
            if extension == ".fasta" or extension == ".FASTA" or extension == ".Fasta":
                self.ui.filepathView_5.setText(self.filepath5[0])
            else:
                choice = QtWidgets.QMessageBox.question(
                    self, 'WARNING!', "Please Choose Fasta File", QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
                if choice == QtWidgets.QMessageBox.Yes:
                    self.browseFile5()
                elif choice == QtWidgets.QMessageBox.No:
                    return
    # Show Result: Pairwise Sequence Alignmen (Global / Local)
    def showResult4(self, filepath1, filepath2):
        if filepath1 != "":
            filepath1 = filepath1[0]
        if filepath2 != "":
            filepath2 = filepath2[0]
        self.ui.msaOutputText_2.clear()
        checkMethod = self.ui.comboBox_4.currentText()
        record1 = SeqIO.parse(filepath1, 'fasta')
        record2 = SeqIO.parse(filepath2, 'fasta')
        Seq1 = ""
        Seq2 = ""
        # Check if user Chose a Global Alignment / Local Alignment methods from a drop down in the GUI 
        if checkMethod == "Global Alignment" or checkMethod == "Local Alignment":
            # If user chose "Global Alignment", then:
            if checkMethod == "Global Alignment":
                for element in record1:
                    Seq1 = str(element.seq)
                for element in record2:
                    Seq2 = str(element.seq)
                # Glocal Alignemnt using globalxx function that takes two sequnces and align them
                globalAlignment = pairwise2.align.globalxx(Seq1, Seq2)
                for index, pair in enumerate(globalAlignment):
                    # Display the Alignment Result on a output Text in GUI
                    self.ui.msaOutputText_2.append(
                        format_alignment((*globalAlignment[index])))

            # If user chose "Local Alignment", then:
            if checkMethod == "Local Alignment":
                for element in record1:
                    Seq1 = str(element.seq)
                for element in record2:
                    Seq2 = str(element.seq)
                # Local Alignemnt using localmx function that takes two sequnces and (match score, mismatch score) and align the two sequences
                localAlignment = pairwise2.align.localmx(Seq1, Seq2, 4, -1)
                for index, pair in enumerate(localAlignment):
                    # Display the Alignment Result on a output Text in GUI
                    self.ui.msaOutputText_2.append(
                        format_alignment((*localAlignment[index])))
        else:
            choice = QtWidgets.QMessageBox.question(
                self, 'WARNING!', "Please Choose Method first", QtWidgets.QMessageBox.Ok)
            if choice == QtWidgets.QMessageBox.Ok:
                return
    # ------------------------------------------------------------------------------------------------------------------ #
    # -------------------------------------------------- End of TAB 1 -------------------------------------------------- #
    # ------------------------------------------------------------------------------------------------------------------ #

def main():
    app = QtWidgets.QApplication(sys.argv)
    application = ApplicationWindow()
    application.show()
    app.exec_()


if __name__ == "__main__":
    main()
