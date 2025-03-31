import subprocess

def run_trimmomatic():
    """
    Trimming low-quality bases using Trimmomatic.
    """
    command = [
        "trimmomatic", "PE", "-phred33", 
        "b1_R1.fastq.gz", "b1_R2.fastq.gz", 
        "b1_R1_trimmed.fastq.gz", "b1_R1_unpaired.fastq.gz", 
        "b1_R2_trimmed.fastq.gz", "b1_R2_unpaired.fastq.gz", 
        "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10", "LEADING:3", "TRAILING:3", 
        "SLIDINGWINDOW:4:15", "MINLEN:36"
    ]
    subprocess.run(command, check=True)

if __name__ == "__main__":
    run_trimmomatic()
