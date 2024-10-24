Proofs of Retrievability with Public Verifiability from Lattices
============================================
These are the experimental evaluation procedures of this paper.

## Description of folders:
* CDH_POR contains the implementation of the publicly verifiable POR scheme based on pairing in [SW13]. The folder 'CDH_POR/src/' contains the implementation codes of CDH-POR scheme. The prerequisites and compilation methods required for testing are shown in the file 'CDH_POR/README'.
* RSA_POR contains the implementation of the publicly verifiable POR scheme based on RSA in [SW13]. The folder 'RSA_POR/src/' contains the implementation codes of RSA-POR scheme. The prerequisites and compilation methods required for testing are shown in the file 'RSA_POR/README'.
* LBPOR contains the implementation of LBPOR scheme. The folder 'LBPOR/src/' contains the implementation codes.The prerequisites and compilation methods required for testing are shown in the file 'LBPOR/README'.
* External_library contains all the external libraries required by the above implementations. 

## External libraries:
* gmp-6.1.2.tar.lz  is the installation package of the GNU multiple precision arithmetic library, which is used to implement various operations of big integers.
* mpfr-4.0.1.tar.gz is the installation package of the GNU MPFR library, which is the prerequisite for installing NFLlib.
* NFLlib-master.zip is the installation package of NTT-based Fast Lattice library, which is used to implement the fast operation of polynomials in ideal lattices.
* openssl-1.1.0l.tar.gz is the installation package of openssl library, in which  implementations of RSAPSS signature and SHA256 are used in our experiments.
* pbc-0.5.14.tar.gz is the installation package of pairing-based cryptography library, which is mainly used to implement various operations on elliptic curves.
