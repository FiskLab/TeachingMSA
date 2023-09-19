//import java.util.Objects;
/**
 * Pair is a class that represents a pair of sequences. It often holds
 * the alignment, the distance, and tree node of the pair of sequences in question
 * @author Nick
 *
 */
public class Pair {
	public Sequence seq1;
	public Sequence seq2;
	public double distanceVal;
	public String[] alignment;
	public String Id;
	double tDist;
	public String tree;
	public String consensus;
	public Pair(Sequence seq1, Sequence seq2){
		this.seq1=seq1;
		this.seq2=seq2;
		this.Id=seq1.Id+"|"+seq2.Id;
		this.alignment=makeAlign(seq1, seq2);
		this.distanceVal=computeDis(alignment);
		
		if(MSA.gimmeTree==true){
			distanceVal=transformDist(seq1, seq2);
		}
		
	}
	/**
	 * updateTree updates the tree in this instance of pair to represent a node in the final
	 * Newick tree. This will be passed to consensus sequence and to other pairs. 
	 */
	public void updateTree(){
		if(seq1.currTree==null&&seq2.currTree==null){
			tree="("+this.seq1.Id+":"+String.format("%.3f",(this.distanceVal/2))+","+
					this.seq2.Id+":"+String.format("%.3f",(this.distanceVal/2))+")";
		}
		else if(seq1.currTree==null&&seq2.currTree!=null){
			tree="("+this.seq1.Id+":"+String.format("%.3f",(this.distanceVal/2))+","+seq2.currTree+
					":"+String.format("%.4f",(this.distanceVal/2))+")";
		}
		else if(seq2.currTree==null&&seq1.currTree!=null){
			tree="("+this.seq2.Id+":"+String.format("%.3f",(this.distanceVal/2))+","+seq1.currTree+
					":"+String.format("%.3f",(this.distanceVal/2))+")";
		}
		else{
			tree="("+this.seq1.currTree+String.format("%.3f",(this.distanceVal/2))+","+seq2.currTree+
					":"+String.format("%.3f",(this.distanceVal/2))+")";
		}
	}
	/**
	 * genConses generates a consensus sequence for a given pair. It picks the
	 * sequence at any given position where one sequence has a gap. If they match, they
	 * keep the sequence, if they mismatch, a random variable is multiplied by the frequency
	 * in all the sequences of that particular residue. The larger value is inserted.
	 * @param alignment
	 * @return
	 */
	public String genConsens(String[] alignment){
			String s1=alignment[0];
			String s2=alignment[1];
			String con="";
			for(int i=0; i<s1.length(); i++){
				if(s2.charAt(i)=='_'&&s1.charAt(i)=='_'){
					s2.replace('_', 'X');
					s1.replace('_', 'X');
				}
				else if(s1.charAt(i)=='_'){
					con+=s2.charAt(i);
				}
				else if(s2.charAt(i)=='_'){
					con+=s1.charAt(i);
				}
				else if(s1.charAt(i)==s2.charAt(i)){
					con+=s1.charAt(i);
				}
				else{
					double val1=MSA.countHash.get(s1.charAt(i)+"")/MSA.countHash.get("total");
					double val2=MSA.countHash.get(s2.charAt(i)+"")/MSA.countHash.get("total");
					val1=val1*Math.random();
					val2=val2*Math.random();
					if(val1>=val2){
						con+=s1.charAt(i);
					}
					else{
						con+=s2.charAt(i);
					}
				}
			}
		
		this.consensus=con;
		return con;
	}
	
	//Left ==[0] in String list
	/**
	 * ComputeDis computes the distance between two sequences(the dissimularity)
	 * 
	 */
	public static double computeDis(String[] alignment){
		String seqAlign1=alignment[0];
		String seqAlign2=alignment[1];
		int totalApply=0;
		int matches=0;
		for(int i=0; i<seqAlign1.length(); i++){
			if(seqAlign1.charAt(i)==seqAlign2.charAt(i)){
				matches+=1;
				totalApply+=1;
			}
			else if(seqAlign1.charAt(i)!='_' && seqAlign2.charAt(i)!='_'){
				totalApply+=1;
			}
		}
		int mismatches=totalApply-matches;
		return (double)mismatches/(double)totalApply;
	}
//	@Override public int hashCode(){
		//return Objects.hash(this.Id);
	//}
	@Override public String toString(){
		return this.Id;
	}
	/**
	 * calls the main MSA to carry out inital alignment calls
	 * @param Seq1
	 * @param Seq2
	 * @return
	 */
	public static String[] makeAlign(Sequence Seq1, Sequence Seq2){
		if(MSA.proteinMode==true){
			return MSA.pairwiseAlignp(Seq1, Seq2);
		}
		else{
			return MSA.pairwiseAlignn(Seq1, Seq2);
		}
	}
	/**
	 * transforms the current distance into a transformed value relative to the outgroup.
	 * @param Seq1
	 * @param Seq2
	 * @return
	 */
	public double transformDist(Sequence Seq1, Sequence Seq2){
		double ds1s2 = 1.1;
		double ds1o = 1.1;
		double ds2o = 1.1;
		double dfinal = 0;
		ds1o = Pair.computeDis(Pair.makeAlign(Seq1, MSA.outgroup));
		ds2o = Pair.computeDis(Pair.makeAlign(Seq2, MSA.outgroup));
		ds1s2 = this.distanceVal;
		dfinal = (((ds1s2 - ds1o - ds2o) / 2.0)+1);

		return dfinal;
	}
}
