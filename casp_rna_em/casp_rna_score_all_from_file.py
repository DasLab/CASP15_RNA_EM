# casp_rna_score_all_from_file
from casp_rna_em.run_metric_programs import score_all_from_file
import argparse


if __name__ == '__main__':
    argParser = argparse.ArgumentParser()
    argParser.add_argument("-f", '--file', type=str, help="file with pdb list")
    argParser.add_argument('--out_file_prefix', type=str, help="")
    argParser.add_argument('--native', type=str, help="")
    argParser.add_argument('--usalign_location', type=str, help="")
    argParser.add_argument('--chimerax_location', type=str, help="")
    argParser.add_argument('--phenix_location', type=str, help="")
    argParser.add_argument('--EM', action='store_true', help="")
    argParser.add_argument('--emmap', type=str, help="")
    argParser.add_argument('--resolution', type=float, help="")
    argParser.add_argument('--threshold', type=float, help="")

    args = argParser.parse_args()

    if args.EM:
        score_all_from_file(file=args.file, out_file_prefix=args.out_file_prefix, native=args.native,
                            usalign_location=args.usalign_location, chimerax_location=args.chimerax_location, EM=args.EM, emmap=args.emmap,
                            resolution=args.resolution, threshold=args.threshold, phenix_location=args.phenix_location)
    else:
        score_all_from_file(file=args.file, out_file_prefix=args.out_file_prefix, native=args.native,
                            usalign_location=args.usalign_location, EM=args.EM, phenix_location=args.phenix_location)
