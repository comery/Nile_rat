#!/usr/bin/env python3
import sys
import os
import re
import argparse

def store_query(faFile):
	dic = {}
	id = ''
	seq = ''
	with open(faFile) as f:
		for line in f:
			line = line.rstrip()
			if line[0] == '>':
				if id != '':
					dic[id[1:]] = len(seq)
				id = line.split()[0]
				seq = ''
			else:
				seq += line
	dic[id[1:]] = len(seq)
	return dic

def parse(inFile, outDir, dic):
	pattern = re.compile(r'(.+?: )(.+?)( :.+?)')
	namePattern = re.compile(r'gene_id (\S+) ; sequence (\S+) ;')

	basename = os.path.basename(inFile)

	gffFile = os.path.join(outDir, '%s.gff' % (basename))
	mutFile = os.path.join(outDir, '%s.mut' % (basename))
	algFile = os.path.join(outDir, '%s.alg' % (basename))

	gff_out = open(gffFile, 'w')
	mut_out = open(mutFile, 'w')
	alg_out = open(algFile, 'w')
	alg_out.write('#GeneID\toriginalLength\tpredict_bg\tpredict_ed\tstrand\tchr\tchr_bg\tchr_ed\talignRate\tscore\tidentity\n')
	
	with open(inFile) as f:
		qryAln = ''
		aln = ''
		refAln = ''
		flag = 0
		a = ''
		same = 0
		
		for line in f:
			line = line.rstrip()

			if line.startswith('  Target range:'):
				f.readline()
				while(1):
					a = f.readline().rstrip()
					if a.startswith('vulgar:'):
						flag = 1
						break
					b = f.readline().rstrip()
					c = f.readline().rstrip()
					f.readline()
					f.readline()
					match = pattern.match(a)
					l = []
					l.append(len(match.group(1)))
					l.append(len(match.group(2)))
					l.append(len(match.group(3)))
					qryAln += a[l[0]:l[0]+l[1]]
					aln += b[l[0]:l[0]+l[1]]
					refAln += c[l[0]:l[0]+l[1]]
					
			if flag == 1: #line.startwith('vulgar:'):
				tmp = a.split()
				qid = tmp[1]
				qbg = int(tmp[2])+1
				qed = int(tmp[3])
				rid = tmp[5]
				rbg = int(tmp[6])
				red = int(tmp[7])
				strand = tmp[8]
				if strand == '-':
					rbg, red = red, rbg
				rbg += 1
				score = int(tmp[9])
				alnRate = (qed-qbg+1)/dic[qid]
				same = aln.count('|')/3
				identity = same/(qed-qbg+1)*100
				alg_out.write('%s\t%i\t%i\t%i\t%s\t%s\t%i\t%i\t%.4f\t%i\t%.2f\n' % (qid, dic[qid], qbg, qed, strand, rid, rbg, red, alnRate, score, identity))
				qryAln = ''
				aln = ''
				refAln = ''
				
				if strand == '+':
					qpos = qbg
					rpos = rbg
				else:
					qpos = qbg
					rpos = red
				for i in range(10, len(tmp), 3):
					typ = tmp[i]
					qlen = int(tmp[i+1])
					rlen = int(tmp[i+2])
					if strand == '+':
						if typ == 'F':
							mut_out.write('%s\t%i\t%i\t%i\t%s\t%s\t%i\t%i\t%.4f\t%i\t%.2f\t%s\t%i\t%i\t%i\t%i\n' % (qid, dic[qid], qbg, qed, strand, rid, rbg, red, alnRate, score, identity, 'F', qpos, qpos, rpos, rpos+rlen-1))
							#mut_out.write('%s\t%i\t%i\t%i\t%s\t%s\t%i\t%i\t$.4f\t%i\t%.2f\t%s\t%i\t%i\t%i\t%i\n' % (qid, dic[qid], qbg, qed, strand, rid, rbg, red, alnRate, score, identity, 'F', qpos, qpos, rpos, rpos+rlen-1))
						qpos += qlen
						rpos += rlen
					else:
						if typ == 'F':
							mut_out.write('%s\t%i\t%i\t%i\t%s\t%s\t%i\t%i\t%.4f\t%i\t%.2f\t%s\t%i\t%i\t%i\t%i\n' % (qid, dic[qid], qbg, qed, strand, rid, rbg, red, alnRate, score, identity, 'F', qpos, qpos, rpos-rlen+1, rpos))
						qpos += qlen
						rpos -= rlen
				flag = 0

			if line == '# seqname source feature start end score strand frame attributes':
				f.readline()
				geneID = ''
				idx = ''
				while(1):
					line = f.readline()
					line = line.rstrip()
					if line == '# --- END OF GFF DUMP ---':
						break
					tmp = line.split('\t')
					if tmp[2] == 'gene':
						tmp[2] = 'mRNA'
						match = namePattern.match(tmp[8])
						idx = match.group(1)
						geneID = match.group(2)
						tmp[8] = 'ID=%s-D%s;' % (geneID, idx)
						out = '\t'.join(tmp)
						gff_out.write(out + '\n')
					elif tmp[2] == 'cds':
						tmp[2] = 'CDS'
						#print(geneID, idx)
						#sys.exit()
						tmp.append('Parent=%s-D%s;' % (geneID, idx))
						out = '\t'.join(tmp)
						gff_out.write(out + '\n')

	gff_out.close()
	mut_out.close()
	alg_out.close()

	

def main():

	parser = argparse.ArgumentParser(description='exonerate output parser')
	parser.add_argument('qry', type=str, help="query.pep")
	parser.add_argument('inFile', type=str, help="out.exonerate")
	parser.add_argument('outDir', type=str, help="outdir")
	
	args = parser.parse_args()

	faFile = args.qry
	inFile = args.inFile
	outDir = args.outDir

	dic = store_query(faFile)

	if not os.path.exists(outDir):
		os.mkdir(outDir)

	parse(inFile, outDir, dic)

if __name__ == '__main__':
	main()

