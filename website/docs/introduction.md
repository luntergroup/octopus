---
id: introduction
title: Introduction
description: Octopus is a haplotype-based variant caller with multiple calling modes.
sidebar_position: 1
---

Octopus is a mapping-based variant caller that implements several calling models within a unified haplotype-aware framework. Each calling model is designed for a particular kind of experimental design. Octopus takes inspiration from particle filtering by constructing a tree of haplotypes and dynamically pruning and extending the tree based on haplotype posterior probabilities in a sequential manner. This allows octopus to implicitly consider all possible haplotypes at a given loci in reasonable time.