import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Scanner;

public class MSA {

	/*
	 * usage is a string printing the usage for the program. Made it private so
	 * it may not be accessed anywhere
	 */
	private static String usage = "Usage: [fasta formatted file]";
	/*
	 * allSeq is an array of all the lines read in from the fasta file. Stores
	 * in order read from file.
	 */
	protected static ArrayList<String> allSeq = new ArrayList<String>();
	/*
	 * formatedSeqList is an array where even numbered indexes are the
	 * description from the fasta file and the even+1 indexes are the
	 * corresponding sequences.
	 */
	protected static ArrayList<String> formatedSeqList = new ArrayList<String>();
	/*
	 * proteinMode stores the mode of the program. If false, the nucleotide
	 * aligner will run. If true, the protein aligner will run. Is accessed by
	 * the function WhichMode().
	 */
	public static boolean proteinMode = false; // False by default
	/*gapPen holds the penalty for introducing a gap*/
	public static int gapPen;
	/*pairwiseAll will hold all the permutations of sequence combinations*/
	public static ArrayList<Pair> pairwiseAll = new ArrayList<Pair>();
	/*gapExtendPen holds the penalty score for an extention of a gap*/
	public static int gapExtendPen;
	/*match holds the value for a match in a nuc alignment*/
	public static int match;
	/*mismatch holds the penalty for a mismatch in a nuc alignment*/
	public static int mismatch;
	/*gimmeTree holds a bool reflecting if the user desires a tree*/
	public static boolean gimmeTree;
	/*pam100 is capable of holding a usable form of any PAM matrix formated like the one
	 * provided*/
	public static String[][] pam100;
	/*seqArray holds all the sequences in the fasta file for use in MSA*/
	public static ArrayList<Sequence> seqArray = new ArrayList<Sequence>();
	/*pamMap makes accessing the pam Matrix faster by storing the index where the 
	 * AA symbols are found in the matrix
	 */
	public static HashMap<String, Integer> pamMap;

	/**
	 * Main takes in a FASTA formated file from the command line.
	 * It splits it up into a form usable by this program and prompts
	 * the user for a series of conditions pertaining to the alignment.
	 * It then calls the appropriate methods to perform to user specs.
	 * @param args
	 * @throws FileNotFoundException
	 */
	public static void main(String[] args) throws FileNotFoundException {
		// TODO Auto-generated method stub
		/*
		 * 
		 * Sequences should be the first and only arguement
		 */
		if (args.length != 1) {// if the number of arguements varies from 1
			System.err.println(usage);// print the usage error
			System.exit(0);// try to exit the program gracefully
		} else {
			// A pretty standard way of reading in from a file.
			BufferedReader buff = new BufferedReader(new FileReader(args[0]));
			String strline;
			try {
				while ((strline = buff.readLine()) != null) {
					allSeq.add(strline);
				}
				buff.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				System.out.println("There was an error with your file. Please"
						+ "check your input files and try again.");
				e.printStackTrace();
			}
			String holdMe = "";
			for (String element : allSeq) {
				if (element.charAt(0) == '>') {
					if (holdMe != "") {
						formatedSeqList.add(holdMe);
						holdMe = "";
					}
					formatedSeqList.add(element);
				}

				else {
					holdMe += element;
				}
			}
			formatedSeqList.add(holdMe);
		}
		formatSeqs();
		generateCounts();

		WhichMode();
		getTree();
		if(gimmeTree==true){
			System.out.println("Using the first Sequence in the FASTA file as the outgroup");
		}
		getGap();
		getExtend();
		if (proteinMode == true) {
			makePAM();
			mkPamHash();
			doMSAp();
		} else {
			getMatch();
			getMismatch();
			doMSAp();
		}
	}
	//outgroupLocation is set permanently to 1.
	public static final int outgroupLocation = 1;
	
	/**
	 * formatSeqs takes the raw input stored in arrays and 
	 * adds just the sequences to seqArray
	 */
	public static void formatSeqs() {
		for (int i = 0; i < formatedSeqList.size(); i++) {
			if (i + 1 == outgroupLocation) {
				outgroup = new Sequence(formatedSeqList.get(i),
						formatedSeqList.get(i + 1), true);
				seqArray.add(new Sequence(formatedSeqList.get(i),
						formatedSeqList.get(i + 1), true));
				i += 1;
			} else {
				seqArray.add(new Sequence(formatedSeqList.get(i),
						formatedSeqList.get(i + 1)));
				i += 1;
			}
		}
	}
	//this is a permanent reference to the Sequence object belonging to the outgroup.
	public static Sequence outgroup;

	/////////Below is optional code to turn selecting the outgroup on for the user
	//////// disabled by default so the user doens't specify an out of bounds outgroup.
	/*public static void getOutgroup() {
		Scanner input = new Scanner(System.in);
		System.out
				.println("Please specify which sequence in your FASTA file is the outgroup: ");
		// input.reset();
		try {
			String response = input.next();
			outgroupLocation = Integer.parseInt(response);
		} catch (NumberFormatException nFE) {
			System.out.println("Integers only, please.");
			input.reset();
			getOutgroup();
		}
		System.out.println("You have selected sequence " + outgroupLocation
				+ " as the outgroup");

	}*/

	/**
	 * getTree prompts the user and asks them if they want a tree.
	 */
	public static void getTree() {
		Scanner scan = new Scanner(System.in);
		System.out.print("UPGMA tree desired? (Newick formatted)(outgroup=1st seq in FASTA)[y,n]");
		char response = scan.next().charAt(0);
		if (response == 'y') {
			gimmeTree = true;
			System.out.println("You got it, boss. One UPGMA tree coming up!");
			// scan.close();
		} else if (response == 'n') {
			gimmeTree = false;
			System.out
					.println("Not feeling a tree today? Maybe next time then...");
			// scan.close();
		} else {
			System.out
					.println("Invalid input: Please try again ('y'==yes, 'n'==no");
			scan.reset();
			getTree();
		}
	}
	/**
	 * makeProMatrix initalizes the matrix to be used for proteins and returns it to 
	 * the user. 
	 * @param Seq1
	 * @param Seq2
	 * @return
	 */
	public static Cell[][] makeProMatrix(Sequence Seq1, Sequence Seq2) {
		int size1 = (Seq1.len + 2);
		int size2 = (Seq2.len + 2);
		Cell[][] pMatrix = new Cell[size1][size2];
		for (int i = 0; i < size1; i++) {
			for (int j = 0; j < size2; j++) {
				if (i == 0) {
					if (j == 0 || j == 1) {
						pMatrix[i][j] = new Cell("+");
					} else {
						pMatrix[i][j] = new Cell(String.valueOf(Seq2.seq
								.charAt(j - 2)));
					}
				} else if (j == 0) {
					if (i == 0 || i == 1) {
						pMatrix[i][j] = new Cell("+");
					} else {
						pMatrix[i][j] = new Cell(String.valueOf(Seq1.seq
								.charAt(i - 2)));
					}
				} else {
					pMatrix[i][j] = new Cell();
				}
			}
		}
		return pMatrix;
	}

	/**
	 * makeNucMatrix initializes the matrix that is filled out during alignment
	 * for the nucleotide alignment. It rejects any non-nucleotide character.
	 * 
	 * @param Seq1
	 * @param Seq2
	 */
	public static Cell[][] makeNucMatrix(Sequence Seq1, Sequence Seq2) {
		int size1 = (Seq1.len + 2);
		int size2 = (Seq2.len + 2);
		Cell[][] nMatrix = new Cell[size1][size2];
		for (int i = 0; i < size1; i++) {
			for (int j = 0; j < size2; j++) {
				if (i == 0) {
					if (j == 0 || j == 1) {
						nMatrix[i][j] = new Cell("+");
					} else {
						if (Seq2.seq.charAt(j - 2) == 'A'
								|| Seq2.seq.charAt(j - 2) == 'T'
								|| Seq2.seq.charAt(j - 2) == 'C'
								|| Seq2.seq.charAt(j - 2) == 'G') {
							nMatrix[i][j] = new Cell(String.valueOf(Seq2.seq
									.charAt(j - 2)));
						} else {
							System.out
									.println("Invalid sequence in your nucleotide file. Please check your file and try again.");
							System.exit(0);
						}
					}
				} else if (j == 0) {
					if (i == 0 || i == 1) {
						nMatrix[i][j] = new Cell("+");
					} else {
						if (Seq1.seq.charAt(i - 2) == 'A'
								|| Seq1.seq.charAt(i - 2) == 'T'
								|| Seq1.seq.charAt(i - 2) == 'C'
								|| Seq1.seq.charAt(i - 2) == 'G') {
							nMatrix[i][j] = new Cell(String.valueOf(Seq1.seq
									.charAt(i - 2)));
						} else {
							System.out
									.println("Invalid sequence in your nucleotide file. Please check your file and try again.");
							System.exit(0);
						}
					}
				} else {
					nMatrix[i][j] = new Cell();
				}
			}

		}
		return nMatrix;
	}
	/**
	 * pairwiseAlignp performs a single pairwise alignment on two input sequences
	 * and returns them as a native array of String.
	 * @param seq1
	 * @param seq2
	 * @return
	 */
	public static String[] pairwiseAlignp(Sequence seq1, Sequence seq2) {
		Cell[][] pMatrix = makeProMatrix(seq1, seq2);
		// printMatrixp(seq1.seq, seq2.seq, pMatrix);
		pMatrix[1][1].selfScore = 0;
		pMatrix[1][1].isStart = true;
		String SeqA = seq1.seq;
		String SeqB = seq2.seq;
		for (int i = 2; i <= SeqA.length() + 1; i++) {
			pMatrix[i][1].selfScore = 0;
			pMatrix[i][1].isGap = true;
			pMatrix[i][1].parent = pMatrix[i - 1][1];
			pMatrix[i][1].pType = "up";
			pMatrix[i][1].selfRow = i;
			pMatrix[i][1].selfCol = 1;

		}
		for (int j = 2; j <= SeqB.length() + 1; j++) {
			pMatrix[1][j].selfScore = 0;
			pMatrix[1][j].isGap = true;
			pMatrix[1][j].parent = pMatrix[1][j - 1];
			pMatrix[1][j].pType = "left";
			pMatrix[1][j].selfRow = 1;
			pMatrix[1][j].selfCol = j;

		}

		pMatrix[2][1].parent = pMatrix[1][1];
		pMatrix[1][2].parent = pMatrix[1][1];

		int diag;
		int up;
		int left;
		int winner;
		pMatrix[1][1].selfScore = 0;
		for (int i = 2; i <= SeqA.length() + 1; i++) {
			for (int j = 2; j <= SeqB.length() + 1; j++) {

				int x = pamMap.get(pMatrix[i][0].seq);
				int y = pamMap.get(pMatrix[0][j].seq);
				int score = Integer.parseInt(pam100[x][y]);
				diag = pMatrix[i - 1][j - 1].selfScore + score;

				if (pMatrix[i - 1][j].isGap == true) {
					up = pMatrix[i - 1][j].selfScore + gapExtendPen;
				} else {
					up = pMatrix[i - 1][j].selfScore + gapPen;
				}
				if (pMatrix[i][j - 1].isGap == true) {
					left = pMatrix[i][j - 1].selfScore + gapExtendPen;
				} else {
					left = pMatrix[i][j - 1].selfScore + gapPen;
				}
				winner = Math.max(Math.max(left, up), diag);
				pMatrix[i][j].selfScore = winner;
				if (diag >= left && diag >= up) {
					pMatrix[i][j].parent = pMatrix[i - 1][j - 1];
					pMatrix[i][j].pType = "diag";
					pMatrix[i - 1][j - 1].child = pMatrix[i][j];
				} else if (left >= diag && left >= up) {
					pMatrix[i][j].parent = pMatrix[i][j - 1];
					pMatrix[i][j].pType = "left";
					pMatrix[i][j - 1].child = pMatrix[i][j];
					pMatrix[i][j].isGap = true;
				} else if (up >= diag && up >= left) {
					pMatrix[i][j].parent = pMatrix[i - 1][j];
					pMatrix[i][j].pType = "up";
					pMatrix[i - 1][j].child = pMatrix[i][j];
					pMatrix[i][j].isGap = true;
				}
				pMatrix[i][j].selfCol = j;
				pMatrix[i][j].selfRow = i;
			}
		}
		return printAlignp(seq1, seq2, pMatrix);

	}
	/**
	 * pairwiseAlignn performs a single alignment on a pair of nucleotide sequences
	 * @param seq1
	 * @param seq2
	 * @return
	 */
	public static String[] pairwiseAlignn(Sequence seq1, Sequence seq2) {
		Cell[][] nMatrix = makeNucMatrix(seq1, seq2);
		nMatrix[1][1].selfScore = 0;
		nMatrix[1][1].isStart = true;
		String SeqA = seq1.seq;
		String SeqB = seq2.seq;
		for (int i = 2; i <= SeqA.length() + 1; i++) {
			nMatrix[i][1].selfScore = 0;
			nMatrix[i][1].isGap = true;
			nMatrix[i][1].parent = nMatrix[i - 1][1];
			nMatrix[i][1].pType = "up";
			nMatrix[i][1].selfRow = i;
			nMatrix[i][1].selfCol = 1;

		}
		for (int j = 2; j <= SeqB.length() + 1; j++) {
			nMatrix[1][j].selfScore = 0;
			nMatrix[1][j].isGap = true;
			nMatrix[1][j].parent = nMatrix[1][j - 1];
			nMatrix[1][j].pType = "left";
			nMatrix[1][j].selfRow = 1;
			nMatrix[1][j].selfCol = j;

		}

		nMatrix[2][1].parent = nMatrix[1][1];
		nMatrix[1][2].parent = nMatrix[1][1];

		int diag;
		int up;
		int left;
		int winner;
		nMatrix[1][1].selfScore = 0;
		for (int i = 2; i <= SeqA.length() + 1; i++) {
			for (int j = 2; j <= SeqB.length() + 1; j++) {
				if (isMatch(nMatrix[i][0], nMatrix[0][j])) {
					diag = nMatrix[i - 1][j - 1].selfScore + match;
				} else {
					diag = nMatrix[i - 1][j - 1].selfScore + mismatch;
				}
				if (nMatrix[i - 1][j].isGap == true) {
					up = nMatrix[i - 1][j].selfScore + gapExtendPen;
				} else {
					up = nMatrix[i - 1][j].selfScore + gapPen;
				}
				if (nMatrix[i][j - 1].isGap == true) {
					left = nMatrix[i][j - 1].selfScore + gapExtendPen;
				} else {
					left = nMatrix[i][j - 1].selfScore + gapPen;
				}
				winner = Math.max(Math.max(left, up), diag);
				nMatrix[i][j].selfScore = winner;
				if (diag >= left && diag >= up) {
					nMatrix[i][j].parent = nMatrix[i - 1][j - 1];
					nMatrix[i][j].pType = "diag";
					nMatrix[i - 1][j - 1].child = nMatrix[i][j];
				} else if (left >= diag && left >= up) {
					nMatrix[i][j].parent = nMatrix[i][j - 1];
					nMatrix[i][j].pType = "left";
					nMatrix[i][j - 1].child = nMatrix[i][j];
					nMatrix[i][j].isGap = true;
				} else if (up >= diag && up >= left) {
					nMatrix[i][j].parent = nMatrix[i - 1][j];
					nMatrix[i][j].pType = "up";
					nMatrix[i - 1][j].child = nMatrix[i][j];
					nMatrix[i][j].isGap = true;
				}
				nMatrix[i][j].selfCol = j;
				nMatrix[i][j].selfRow = i;
			}
		}
		return printAlignp(seq1, seq2, nMatrix);

	}
	/**
	 * prints the constructed matrix. Used for debugging, mostly.
	 * @param Seq1
	 * @param Seq2
	 * @param pMatrix
	 */
	public static void printMatrixp(String Seq1, String Seq2, Cell[][] pMatrix) {
		int size1 = (Seq1.length() + 2);
		int size2 = (Seq2.length() + 2);
		for (int i = 0; i < size1; i++) {
			for (int j = 0; j < size2; j++) {

				System.out.print(pMatrix[i][j]);
				System.out.print("\t");
			}
			System.out.println();
		}
	}
	/**
	 * printAlignp is called by the pairwise aligner to pick the max sequence
	 * (and call breakties to break ties), trace back the matrix, and correct for
	 * any gap-related adjustments.
	 * @param Seq1
	 * @param Seq2
	 * @param pMatrix
	 * @return
	 */
	public static String[] printAlignp(Sequence Seq1, Sequence Seq2,
			Cell[][] pMatrix) {
		String AlignLeft = "";
		String AlignTop = "";
		ArrayList<Cell> max = new ArrayList<Cell>();
		for (int i = 1; i <= Seq2.len + 1; i++) {
			if (i == 1) {
				max.add(pMatrix[Seq1.len + 1][i]);
			}
			if (pMatrix[Seq1.len + 1][i].selfScore > max.get(0).selfScore) {
				max.clear();
				max.add(pMatrix[Seq1.len + 1][i]);
			}
			if (pMatrix[Seq1.len + 1][i].selfScore == max.get(0).selfScore) {
				max.add(pMatrix[Seq1.len + 1][i]);
			}
		}
		for (int j = 1; j <= Seq1.len + 1; j++) {
			if (pMatrix[j][Seq2.len + 1].selfScore > max.get(0).selfScore) {
				max.clear();
				max.add(pMatrix[j][Seq2.len + 1]);
			}
			if (pMatrix[j][Seq2.len + 1].selfScore == max.get(0).selfScore) {
				max.add(pMatrix[j][Seq2.len + 1]);
			}
		}
		String AlignTopCorrection = "";
		String AlignLeftCorrection = "";
		ArrayList<String> winnersLeft = new ArrayList<String>();
		ArrayList<String> winnersTop = new ArrayList<String>();
		int i = 0;
		for (Cell align : max) {
			AlignTop = "";
			AlignLeft = "";
			AlignTopCorrection = "";
			AlignLeftCorrection = "";
			int correctRow = pMatrix.length;
			int currRow = align.selfRow;
			int correctCol = pMatrix[1].length;
			int currCol = align.selfCol;
			while (currRow < correctRow) {
				AlignTopCorrection = AlignTopCorrection + '_';
				AlignLeftCorrection = AlignLeftCorrection
						+ Seq1.seq.charAt(currRow - 2);
				currRow++;
			}
			while (currCol < correctCol) {
				AlignLeftCorrection = AlignLeftCorrection + "_";
				AlignTopCorrection = AlignTopCorrection
						+ Seq2.seq.charAt(currCol - 2);
				currCol++;
			}
			Cell curr = align;
			while (curr.isStart != true) {
				if (curr.pType == "up") {
					AlignTop = "_" + AlignTop;
					AlignLeft = pMatrix[curr.selfRow][0].seq + AlignLeft;

				} else if (curr.pType == "left") {
					AlignLeft = "_" + AlignLeft;
					AlignTop = pMatrix[0][curr.selfCol].seq + AlignTop;
				} else if (curr.pType == "diag") {

					AlignLeft = pMatrix[curr.selfRow][0].seq + AlignLeft;
					AlignTop = pMatrix[0][curr.selfCol].seq + AlignTop;// }
				}
				curr = curr.parent;
			}

			if (AlignLeftCorrection.length() - 1 >= 1) {
				AlignLeft += AlignLeftCorrection.substring(1,
						AlignLeftCorrection.length() - 1);
			}
			if (AlignTopCorrection.length() - 1 >= 1) {
				AlignTop += AlignTopCorrection.substring(1,
						AlignTopCorrection.length() - 1);
			}
			winnersLeft.add(AlignLeft);
			winnersTop.add(AlignTop);
		}
		String[] finalWinners = breakties(winnersLeft, winnersTop);
		AlignTop = finalWinners[1];
		AlignLeft = finalWinners[0];
		//String printtwo = formatedSeqList.get(2) + "          ";
		//String printone = formatedSeqList.get(0) + "          ";
		/*
		 * int lineStart=0; int lineEnd=0; while (lineEnd<AlignTop.length()){
		 * lineEnd+=50; if(lineEnd>AlignTop.length()){
		 * lineEnd=AlignTop.length(); } System.out.print(printtwo.substring(0,
		 * 11) + ": "); System.out.println(AlignTop.substring(lineStart,
		 * lineEnd));// Seq 1 is on left System.out.print("             ");
		 * while (i < lineEnd) { System.out.print("|"); i++; }
		 * System.out.println(); System.out.print(printone.substring(0, 11) +
		 * ": "); System.out.println(AlignLeft.substring(lineStart,lineEnd));//
		 * Seq 2 is on top lineStart=lineEnd; }
		 */
		// }
		String[] temp = new String[2];
		temp[0] = AlignLeft;
		temp[1] = AlignTop;
		return temp;
	}
	/**
	 * breakties breaks ties in maximum matrix score by comparing the number of internal gaps
	 * and picking the one with the least.
	 */
	public static String[] breakties(ArrayList<String> winnersLeft,
			ArrayList<String> winnersTop) {
		ArrayList<String> winnersBoth = new ArrayList<String>();
		for (int z = 0; z < winnersLeft.size(); z++) {
			winnersBoth.add(winnersLeft.get(z) + winnersTop.get(z));
		}
		boolean frontgap = true;
		boolean endgap = false;
		int frontgapcount = 0;
		int endgapcount = 0;
		int allgapcount = 0;
		String currwinner = "";
		double currwinnerscore = 0;
		for (String s : winnersBoth) {
			for (int i = 0; i < s.length(); i++) {
				if (s.charAt(i) == '_') {
					allgapcount += 1;
					if (frontgap == true) {
						frontgapcount += 1;
					}
				} else {
					frontgap = false;
				}
			}
			for (int j = s.length() - 1; j >= 0; j--) {
				if (s.charAt(j) == '_') {
					endgapcount += 1;
				} else {
					break;
				}
			}
			if (currwinner == "") {
				currwinner = s;
				currwinnerscore = ((allgapcount - frontgapcount - endgapcount) / (double) allgapcount);
			} else if (currwinnerscore > ((allgapcount - frontgapcount - endgapcount) / (double) allgapcount)) {
				currwinnerscore = ((allgapcount - frontgapcount - endgapcount) / (double) allgapcount);
				currwinner = s;
			}
		}
		int indexofwinners;
		indexofwinners = winnersBoth.indexOf(currwinner);
		String[] finalWinners = new String[2];
		finalWinners[0] = winnersLeft.get(indexofwinners);
		finalWinners[1] = winnersTop.get(indexofwinners);
		return finalWinners;
	}
	/**
	 * Make initial alignments makes all of the inital alignments (in the form 
	 * of Pair Objects) and returns the one with the lowest disimularity.
	 * @return
	 */
	public static Pair mkInitAlignments() {
		double smallestDist = Double.POSITIVE_INFINITY;
		Pair currentWin = null;
			for (int i = 0; i < seqArray.size(); i++) {
				for (int j = i + 1; j < seqArray.size(); j++) {
					Pair temp = new Pair(seqArray.get(i), seqArray.get(j));
					pairwiseAll.add(temp);
					if (temp.distanceVal < smallestDist) {
						currentWin = temp;
						smallestDist = temp.distanceVal;
					}

				}
			}
		return currentWin;
	}
	/**
	 * alignM performs all of the serial pairwise alignments until only one sequence remains
	 * @param first
	 * @param tempSeqList
	 * @param tempPairList
	 * @return
	 */
	public static Sequence alignM(Pair first, ArrayList<Sequence> tempSeqList,
			ArrayList<Pair> tempPairList) {
		if(gimmeTree==true){
			first.updateTree();
		}
		Pair tempPair = null;
		Sequence conSeq = null;
		while (tempSeqList.size() > 1) {
			first.genConsens(first.alignment);
			tempPairList.remove(first);
			Iterator<Pair> piter = tempPairList.iterator();
			while (piter.hasNext()) {
				Pair curr = piter.next();
				if (curr.seq1 == first.seq1 || curr.seq1 == first.seq2) {
					piter.remove();
				}
			}
			Iterator<Pair> piterDos = tempPairList.iterator();
			while (piterDos.hasNext()) {
				Pair curr = piterDos.next();
				if (curr.seq2 == first.seq1 || curr.seq2 == first.seq2) {
					piterDos.remove();
				}
			}
			Pair distWin = null;
			double bestDist = Double.POSITIVE_INFINITY;
			tempSeqList.remove(first.seq1);
			tempSeqList.remove(first.seq2);
			conSeq = new Sequence(first.Id, first.consensus, first.tree);
			for (int i = 0; i < tempSeqList.size(); i++) {
				tempPair = new Pair(tempSeqList.get(i), conSeq);
				tempPairList.add(tempPair);
			}
			for (int q = 0; q < tempPairList.size(); q++) {
				if (tempPairList.get(q).distanceVal < bestDist) {
					bestDist = tempPair.distanceVal;
					distWin = tempPair;
				}
			}
			tempSeqList.add(conSeq);
			if(gimmeTree==true){
				if(tempSeqList.size()>1){
				distWin.updateTree();
				}
			}
			first = distWin;
			

		}
		return tempSeqList.get(0);
	}
	
	 //get rid is a debugging array storing and fixing troublesome indexes 
	public static ArrayList<Integer> getRid = new ArrayList<Integer>();
	/**
	 * finalConsensus takes the final alignments to the consensus and generates 
	 * the final consensus to be presented
	 * @param finalAlignments
	 * @return
	 */
	public static String finalConsensus(String[] finalAlignments) {
		String current = "";
		int count_ = 0;
		// System.out.println(finalAlignments[0]+'\n'+finalAlignments[1]);
		int MaxLen = 0;
		for (int q = 1; q < finalAlignments.length; q += 2) {
			if (finalAlignments[q].length() >= MaxLen) {
				MaxLen = finalAlignments[q].length();
			}
		}
		for (int i = 1; i < finalAlignments.length; i += 2) {
			while (finalAlignments[i].length() < MaxLen) {
				finalAlignments[i] += "_";
			}
		}
		for (int j = 0; j < finalAlignments[1].length(); j++) {
			boolean same = true;
			count_ = 0;
			HashMap<String, Integer> tempHash = new HashMap<String, Integer>();
			for (int i = 1; i < finalAlignments.length; i += 2) {
				if (finalAlignments[i].charAt(j) == '_') {
					same = false;
					count_ += 1;
					if (count_ == finalAlignments.length / 2) {
						count_ = 0;
						getRid.add(j);
					}
				}
				if (finalAlignments[i].charAt(j) != '_') {
					count_ = 0;
					if (!tempHash
							.containsKey(finalAlignments[i].charAt(j) + "")) {
						tempHash.put(finalAlignments[i].charAt(j) + "", 1);
						if (tempHash.size() >= 2) {
							same = false;
						}
					} else {
						tempHash.put(
								finalAlignments[i].charAt(j) + "",
								tempHash.get(finalAlignments[i].charAt(j) + "") + 1);
					}
				}
			}
			if (same == true) {
				current = current + "*";
			} else if (tempHash.size() > 1) {
				int maxValInMap = (Collections.max(tempHash.values()));
				ArrayList<String> temp = new ArrayList<String>();
				for (String entry : tempHash.keySet()) {
					if (tempHash.get(entry) == maxValInMap) {
						temp.add(entry);
					}
				}
				if (tempHash.size() > 1 && temp.size() == 1) {
					current = current
							+ (Character.toString(Character.toLowerCase((temp
									.get(0).charAt(0)))));
				} else if (tempHash.size() > 1 && temp.size() > 1) {
					current = current
							+ (Character.toString(Character.toLowerCase((temp
									.get(1).charAt(0)))));
				}
			} else if (tempHash.size() == 1 && same == false) {
				current = current + ("" + tempHash.keySet()).substring(1, 2);
			}

		}
		return current;
	}
	/**
	 * doMSAp calls mkInitAlignments, alignM, and finalConsensus to perform the bulk of 
	 * the alignment process.It prints the alignments, consensus, and, if applicable, tree.
	 */
	public static void doMSAp() {
		Sequence conSeq = null;
		Pair first = mkInitAlignments();
		ArrayList<Sequence> tempSeqList = new ArrayList<Sequence>(seqArray);
		ArrayList<Pair> tempPairList = new ArrayList<Pair>(pairwiseAll);
		conSeq = alignM(first, tempSeqList, tempPairList);
		conSeq.Id = "Consensus";
		String[] temp = new String[2];
		String[] finalAlign = new String[(2 * (seqArray.size()))];
		for (int i = 0; i < seqArray.size(); i++) {
			temp = Pair.makeAlign(conSeq, seqArray.get(i));
			finalAlign[2 * i] = seqArray.get(i).Id;
			finalAlign[2 * i + 1] = temp[1];
		}
		String printConsens = finalConsensus(finalAlign);
		for (int j = 0; j < finalAlign.length; j += 2) {
			ArrayList<Integer> tempL = new ArrayList<Integer>(getRid);
			//System.out.println(getRid);
			while (!tempL.isEmpty()) {
				finalAlign[j + 1] = finalAlign[j + 1].substring(0,
						tempL.get(tempL.size() - 1))
						+ finalAlign[j + 1]
								.substring(tempL.get(tempL.size() - 1) + 1);
				tempL.remove(tempL.size() - 1);
			}
			System.out.print((finalAlign[j] + "          ").substring(0, 11));
			System.out.println(": " + finalAlign[j + 1]);
		}
		System.out.print("Consensus  ");
		System.out.println(": " + printConsens);
		if(gimmeTree==true){
			System.out.println('\n'+"Here is the Newick Formated UPGMA Tree:");
			System.out.println(conSeq.currTree);
		}
	}
	/**
	 * makes the hashmap for the pam matrix
	 */
	public static void mkPamHash() {
		pamMap = new HashMap<String, Integer>();
		for (int i = 0; i < 21; i++) {
			pamMap.put(pam100[0][i], i);
		}

	}

	public static String getPamFileName() {
		Scanner input = new Scanner(System.in);
		System.out.println("Pam Matrix file?: ");
		// input.reset();
		if (input.hasNext()) {
			String response = input.next();
			return response;
		} else {
			System.out.println("Please enter a file");
			input.reset();
			getPamFileName();
		}
		return "No file";
	}

	public static void makePAM() throws FileNotFoundException {
		String file = getPamFileName();
		Scanner scan = new Scanner(new BufferedReader(new FileReader(file)));
		String[] temp = scan.nextLine().trim().split("\\s+");
		pam100 = new String[21][21];
		pam100[0][0] = "-";
		// System.out.print(pam100[0][0]+" ");
		for (int q = 1; q < 21; q++) {
			pam100[0][q] = temp[q - 1];
			// System.out.print(pam100[0][q]+" ");
		}
		for (int i = 1; i < 21; i++) {
			// System.out.println();
			temp = scan.nextLine().trim().split("\\s+");
			// System.out.println(temp);
			for (int j = 0; j < 21; j++) {
				pam100[i][j] = temp[j];
				// System.out.print(pam100[i][j]+ " ");
			}
		}

	}

	public static void WhichMode() {
		Scanner scan = new Scanner(System.in);
		System.out
				.print("Enter 'n' for nucleotide or 'p' for protein alignment: ");
		char response = scan.next().charAt(0);
		if (response == 'p') {
			proteinMode = true;
			System.out.println("You have selected protein alignment...");
			// scan.close();
		} else if (response == 'n') {
			proteinMode = false;
			System.out.println("You have selected nucleotide alignment...");
			// scan.close();
		} else {
			System.out.println("Invalid input: Please try again");
			scan.reset();
			WhichMode();
		}
	}

	public static void getGap() {
		Scanner input = new Scanner(System.in);
		System.out.println("Gap Penalty Score?: ");
		// input.reset();
		try {
			String response = input.next();
			gapPen = Integer.parseInt(response);
		} catch (NumberFormatException nFE) {
			System.out.println("Integers only, please.");
			input.reset();
			getGap();
		}
		System.out.println("You have selected a gap penalty of " + gapPen);
	}

	public static void getExtend() {
		Scanner input = new Scanner(System.in);
		System.out.println("Gap Extention Penalty Score?: ");
		// input.reset();
		try {
			String response = input.next();
			gapExtendPen = Integer.parseInt(response);
		} catch (NumberFormatException nFE) {
			System.out.println("Integers only, please.");
			input.reset();
			getExtend();
		}
		System.out.println("You have selected a gap extention penalty of "
				+ gapExtendPen);
	}

	public static void getMatch() {
		Scanner input = new Scanner(System.in);
		System.out.println("Match Score?: ");
		// input.reset();
		try {
			String response = input.next();
			match = Integer.parseInt(response);
		} catch (NumberFormatException nFE) {
			System.out.println("Integers only, please.");
			input.reset();
			getMatch();
		}
		System.out.println("You have selected a Match Score of " + match);
	}

	/**
	 * get misMatch prompts the user for the mismatch penalty
	 */
	public static void getMismatch() {
		Scanner input = new Scanner(System.in);
		System.out.println("Mismatch Penalty Score?: ");
		// input.reset();
		try {
			String response = input.next();
			mismatch = Integer.parseInt(response);
		} catch (NumberFormatException nFE) {
			System.out.println("Integers only, please.");
			input.reset();
			getMismatch();
		}
		System.out.println("You have selected a mismatch penalty of "
				+ mismatch);
	}

	public static boolean isMatch(Cell left, Cell top) {
		if (!left.seq.equals(top.seq)) {
			return false;
		}
		return true;
	}
	/**
	 * generate counts counts the amounts of different building blocks and stores them
	 * for use in breaking ties in consensus.
	 */
	public static HashMap<String, Integer> countHash = new HashMap<String, Integer>();

	public static void generateCounts() {
		int total = 0;
		for (int i = 1; i < formatedSeqList.size(); i += 2) {
			String tempSeq = formatedSeqList.get(i);
			for (int j = 0; j < tempSeq.length(); j++) {
				if (!countHash.containsKey(tempSeq.charAt(j) + "")) {
					countHash.put(tempSeq.charAt(j) + "", 1);
					total += 1;
				} else {
					countHash.put(tempSeq.charAt(j) + "",
							countHash.get(tempSeq.charAt(j) + "") + 1);
					total += 1;
				}
			}
		}
		countHash.put("total", total);
	}

}