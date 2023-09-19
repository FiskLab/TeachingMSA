
/**
 * Sequence is an object representing a single sequence. It can also
 * hold the information about the tree it represents, if it is a consensus
 * sequence. It has 3 different constructors! 
 * @author Nick
 *
 */
public class Sequence {
	public String Id;
	public String seq;
	public String currTree;
	public boolean isOutgroup;
	public int len;
	public Sequence(String Id, String seq){
		this.Id=Id;
		this.seq=seq;
		this.len=seq.length();
		this.isOutgroup=false;
	}
	public Sequence(String Id, String seq, String currTree){
		this.Id=Id;
		this.seq=seq;
		this.len=seq.length();
		this.isOutgroup=false;
		this.currTree=currTree;
	}
	public Sequence(String Id, String seq, Boolean holdmyBoolson){
		this.Id=Id;
		this.seq=seq;
		this.len=seq.length();
		this.isOutgroup=holdmyBoolson;

	}
	@Override public String toString(){
		return this.seq;
	}
}
