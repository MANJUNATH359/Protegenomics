from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable

INPUT_FASTA = "AmoebaDB-68_EhistolyticaHM1IMSS_Genome.fasta"
OUTPUT_FASTA = "Final_ORFs.fasta"
OUTPUT_LOG = "Final_ORFs.log"

MIN_AA = 10
START_CODON = "ATG"

table = CodonTable.unambiguous_dna_by_name["Standard"]
genetic_code = table.forward_table
stop_codons = set(table.stop_codons)


def scan_frame(seq, offset, strand, accession, seq_len, fa, logf):
    pos = offset
    found_first_orf = False
    orf_id = 0

    while pos + 3 <= seq_len:
        codon = seq[pos:pos+3]

        # ---- FIRST ORF MUST START WITH ATG ----
        if not found_first_orf:
            if codon != START_CODON:
                pos += 3
                continue
            aa = ["M"]
            start_nt = pos + 1
            cur = pos + 3

        # ---- UNBIASED ORFs ----
        else:
            if "N" in codon or codon in stop_codons:
                pos += 3
                continue
            aa = [genetic_code.get(codon, "X")]
            start_nt = pos + 1
            cur = pos + 3

        end_nt = None
        stop_reason = "END"

        while cur + 3 <= seq_len:
            nxt = seq[cur:cur+3]

            if "N" in nxt:
                stop_reason = "N"
                end_nt = cur - 3
                break

            if nxt in stop_codons:
                stop_reason = nxt
                end_nt = cur
                cur += 3
                break

            aa.append(genetic_code.get(nxt, "X"))
            cur += 3

        if end_nt is None:
            end_nt = cur - 3

        if len(aa) >= MIN_AA:
            orf_id += 1

            if strand == "+":
                g_start = start_nt
                g_end = end_nt + 3
                frame_label = f"+{offset+1}"
            else:
                g_start = seq_len - (end_nt + 3) + 1
                g_end = seq_len - (start_nt - 1)
                frame_label = f"-{offset+1}"

            header = f"{accession}_{frame_label}_{g_start}-{g_end}_ORF{orf_id}"

            fa.write(f">{header}\n{''.join(aa)}\n")
            logf.write(
                f"{header}\t{accession}\t{frame_label}\t"
                f"{g_start}\t{g_end}\t{stop_reason}\t{len(aa)}\n"
            )

            found_first_orf = True
            pos = cur
        else:
            pos += 3


def run():
    with open(OUTPUT_FASTA, "w") as fa, open(OUTPUT_LOG, "w") as logf:
        logf.write("ORF_ID\tACCESSION\tFRAME\tSTART\tEND\tSTOP\tAA_LEN\n")

        for rec in SeqIO.parse(INPUT_FASTA, "fasta"):
            acc = rec.id
            seq = str(rec.seq).upper()
            rc = str(Seq(seq).reverse_complement())
            L = len(seq)

            # +1 +2 +3
            for offset in (0, 1, 2):
                scan_frame(seq, offset, "+", acc, L, fa, logf)

            # -1 -2 -3
            for offset in (0, 1, 2):
                scan_frame(rc, offset, "-", acc, L, fa, logf)

    print("Six-frame ORF translation completed correctly")


if __name__ == "__main__":
    run()