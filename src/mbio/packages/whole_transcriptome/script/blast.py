'''''
@author: JiuYu
'''


def score(a, b):  # scoring function
    score = 0
    lst = ['AC', 'GT', 'CA', 'TG']
    if a == b:
        score += 2
    elif a + b in lst:
        score += -5
    else:
        score += -7
    return score


def BLAST(seq1, seq2):  # Basic Local Alignment Search Tool
    l1 = len(seq1)
    l2 = len(seq2)
    GAP = -5  # -5 for any gap
    scores = []
    point = []

    for j in range(l2 + 1):
        if j == 0:
            line1 = [0]
            line2 = [0]
            for i in range(1, l1 + 1):
                line1.append(GAP * i)
                line2.append(2)
        else:
            line1 = []
            line2 = []
            line1.append(GAP * j)
            line2.append(3)
        scores.append(line1)
        point.append(line2)

    # fill the blank of scores and point
    for j in range(1, l2 + 1):
        letter2 = seq2[j - 1]
        for i in range(1, l1 + 1):
            letter1 = seq1[i - 1]
            diagonal_score = score(letter1, letter2) + scores[j - 1][i - 1]
            left_score = GAP + scores[j][i - 1]
            up_score = GAP + scores[j - 1][i]
            max_score = max(diagonal_score, left_score, up_score)
            scores[j].append(max_score)

            if scores[j][i] == diagonal_score:
                point[j].append(1)
            elif scores[j][i] == left_score:
                point[j].append(2)
            else:
                point[j].append(3)

    # trace back
    alignment1 = ''
    alignment2 = ''
    i = l2
    j = l1
    print 'scores =', scores[i][j]
    while True:
        if point[i][j] == 0:
            break
        elif point[i][j] == 1:
            alignment1 += seq1[j - 1]
            alignment2 += seq2[i - 1]
            i -= 1
            j -= 1
        elif point[i][j] == 2:
            alignment1 += seq1[j - 1]
            alignment2 += '-'
            j -= 1
        else:
            alignment1 += '-'
            alignment2 += seq2[i - 1]
            i -= 1

    # reverse alignment
    alignment1 = alignment1[::-1]
    alignment2 = alignment2[::-1]
    print 'The best alignment:'
    print alignment1
    print alignment2


seq1 = raw_input('Please input your first sequences:\n')
seq2 = raw_input('input second sequences:\n')
BLAST(seq1, seq2)
