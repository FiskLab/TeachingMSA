
/**
 * Cell represents an individual unit of a scoring matrix. It is mostly a container
 * for data. 
 * @author Nick
 *
 */
public class Cell {
	public String pType;
	public int selfRow;
	public int selfCol;
	public int selfScore;
	public boolean isSeq;
	public String seq;
	public Cell parent;
	public boolean isGap;
	public Cell child;
	public boolean isStart;
	
	public Cell(){
		this.isSeq=false;
		//this.isResolved=false;
		//this.timesVisited=0;
		this.isGap=false;
		this.isStart=false;
	}
	
	public Cell(String seq){
		//this.isResolved=false;
		//this.timesVisited=0;
		this.seq=seq;
		this.isSeq=true;
		this.isStart=false;
	}
	@Override public String toString(){
		if(this.isSeq==true){
			return this.seq;
		}
		else{
			return String.valueOf(this.selfScore);
		}
	}
	
}
