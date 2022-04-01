# -*- coding: utf-8 -*-


class ListNode(object):
    def __init__(self, x):
        self.val = x
        self.next = None


if __name__ == '__main__':
    head = ListNode(None)
    n1 = ListNode(1)
    n2 = ListNode(2)
    n3 = ListNode(3)
    head.next = n1
    n1.next = n2
    n2.next = n3
    print(head)
