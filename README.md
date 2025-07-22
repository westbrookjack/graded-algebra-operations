# Graded Algebra Operations in Macaulay2

This repository contains Macaulay2 scripts for computing presentations of:

- **Segre products** of graded rings  
- **Veronese subrings** of a given graded ring

These tools are intended for use in commutative algebra and algebraic geometry research, particularly in the study of Hilbert–Kunz multiplicities and projective coordinate rings.

## Contents

- `SegreProduct.m2` —  
  Computes a presentation for the Segre product of two finitely generated standard graded \( k \)-algebras.  
  The resulting ring is given with generators and relations, suitable for further exploration (e.g. computing Hilbert series, multiplicities, or invariants).

- `VeroneseSubring.m2` —  
  Constructs the \( d \)-th Veronese subring of a standard graded \( k \)-algebra by producing a presentation in new generators of degree 1.  
  Compatible with downstream computations like syzygies or resolutions.

## Usage

To use these scripts, load them into a Macaulay2 session:

```macaulay2
load "SegreProduct.m2"
load "VeroneseSubring.m2"
```
## Background & Motivation

These tools were developed as part of an REU project at the University of Michigan-Ann Arbor on Hilbert-Kunz multiplicities of Segre products under the guidance of Dr. Austyn Simpson in 2024.

##License

No formal license is currently specified. If you'd like to use or adapt these scripts in your own work, feel free to reach out.

##Contact

Jack Westbrook
