import random
import typing

import pytest


def make_ansatz(s: str, n: int, i: int) -> typing.Optional[str]:
    if not 0 <= i + n <= len(s):
        return None
    else:
        return s[: i + n] + s[i:]


def score_alignments(first: str, second: str) -> typing.Iterable[int]:
    n = len(second) - len(first)
    for i in range(len(first)):
        ansatz = make_ansatz(first, n, i)
        if ansatz is not None:
            indices = [*range(len(first))]
            yield (
                sum(a == b for a, b in zip(ansatz, second)),  # score
                make_ansatz(indices, n, i),  # alignment
            )


def find_best_alignment(first: str, second: str) -> typing.Optional[str]:
    scored_alignments = [*score_alignments(first, second)]
    best_score = max(score for score, _alignment in scored_alignments)
    rng = random.Random(first + second)
    best_scoring = [
        alignment
        for score, alignment in scored_alignments
        if score == best_score
    ]
    return rng.choice(best_scoring)


def test_find_best_alignment_basic1():
    first = "abcdef"
    second = "abcbcdef"
    result = find_best_alignment(first, second)
    assert result == [0, 1, 2, 1, 2, 3, 4, 5]


def test_find_best_alignment_basic2():
    first = "abcdef"
    second = "adef"
    result = find_best_alignment(first, second)
    assert result == [0, 3, 4, 5]

def test_find_best_alignment_basic3():
    first = "abcdef"
    second = "abydef"
    result = find_best_alignment(first, second)
    assert result == [0, 1, 2, 3, 4, 5]
