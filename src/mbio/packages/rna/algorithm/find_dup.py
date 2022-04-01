# -*- coding: utf-8 -*


def find_dup(arr):
    i = 0
    while True:
        if i == len(arr):
            return -1
        if arr[i] == i:
            i += 1
            continue
        v = arr[i]
        if arr[i] == arr[v]:
            break
        arr[i], arr[v] = arr[v], arr[i]
    return arr[i]


if __name__ == '__main__':
    array = [3, 2, 4, 1, 3, 5]
    dup = find_dup(array)
    print(dup)
