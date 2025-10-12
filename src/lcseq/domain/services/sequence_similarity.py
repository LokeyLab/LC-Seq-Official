"""Sequence similarity analysis domain service.

This module provides domain services for computing sequence similarity metrics,
including edit distance calculations between sequences.

Domain Responsibility:
    - Compute similarity/distance metrics between sequences
    - Provide various distance algorithms (Levenshtein, etc.)
    - Pure algorithmic logic with no presentation concerns

Not Responsible For:
    - Sorting sequences for display (that's presentation)
    - Visualization or plotting
    - File I/O or data loading
"""

from typing import List


class SequenceSimilarityAnalyzer:
    """Domain service for sequence similarity computation.

    This service provides algorithms for computing similarity and distance
    metrics between text sequences. All methods are deterministic and side-effect free.

    Examples:
        >>> analyzer = SequenceSimilarityAnalyzer()
        >>> distance = analyzer.levenshtein_distance("kitten", "sitting")
        >>> distance
        3
        >>> distance = analyzer.levenshtein_distance("ABC-DEF-GHI", "ABC-DEF-XYZ")
        >>> distance
        3
    """

    @staticmethod
    def levenshtein_distance(s1: str, s2: str) -> int:
        """Compute Levenshtein (edit) distance between two strings.

        The Levenshtein distance is the minimum number of single-character edits
        (insertions, deletions, or substitutions) required to transform one string
        into another. This is also known as edit distance.

        Algorithm:
            Uses dynamic programming with O(min(m,n)) space complexity.
            Time complexity is O(m*n) where m and n are string lengths.

        Args:
            s1: First string
            s2: Second string

        Returns:
            Non-negative integer representing the minimum edit distance

        Examples:
            >>> SequenceSimilarityAnalyzer.levenshtein_distance("", "")
            0
            >>> SequenceSimilarityAnalyzer.levenshtein_distance("", "abc")
            3
            >>> SequenceSimilarityAnalyzer.levenshtein_distance("kitten", "sitting")
            3
            >>> SequenceSimilarityAnalyzer.levenshtein_distance("Saturday", "Sunday")
            3

        References:
            Levenshtein, Vladimir I. (1966). "Binary codes capable of correcting
            deletions, insertions, and reversals". Soviet Physics Doklady. 10 (8): 707–710.
        """
        # Optimization: ensure s1 is the longer string to minimize space
        if len(s1) < len(s2):
            return SequenceSimilarityAnalyzer.levenshtein_distance(s2, s1)

        # Base case: if s2 is empty, distance is length of s1
        if len(s2) == 0:
            return len(s1)

        # Dynamic programming: only need previous row, not full matrix
        previous_row = list(range(len(s2) + 1))

        for i, c1 in enumerate(s1):
            current_row = [i + 1]  # First element is always i+1

            for j, c2 in enumerate(s2):
                # Cost of operations:
                insertions = previous_row[j + 1] + 1       # Insert c1
                deletions = current_row[j] + 1             # Delete c2
                substitutions = previous_row[j] + (c1 != c2)  # Substitute if different

                current_row.append(min(insertions, deletions, substitutions))

            previous_row = current_row

        return previous_row[-1]

    def compute_pairwise_distances(self, sequences: List[str]) -> List[List[int]]:
        """Compute pairwise Levenshtein distances between all sequences.

        Args:
            sequences: List of strings to compare

        Returns:
            Square matrix (as list of lists) where matrix[i][j] is the
            distance between sequences[i] and sequences[j]

        Examples:
            >>> analyzer = SequenceSimilarityAnalyzer()
            >>> seqs = ["ABC", "ABD", "XYZ"]
            >>> matrix = analyzer.compute_pairwise_distances(seqs)
            >>> matrix[0][1]  # Distance between "ABC" and "ABD"
            1
            >>> matrix[0][2]  # Distance between "ABC" and "XYZ"
            3
        """
        n = len(sequences)
        distance_matrix = [[0] * n for _ in range(n)]

        for i in range(n):
            for j in range(i + 1, n):
                dist = self.levenshtein_distance(sequences[i], sequences[j])
                distance_matrix[i][j] = dist
                distance_matrix[j][i] = dist

        return distance_matrix
