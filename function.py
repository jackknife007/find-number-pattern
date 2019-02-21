# -*- coding: utf-8 -*-

import numpy as np
from functools import reduce


class Sequence(object):

    def __init__(self, operator=None, sub=None, status=None, isprint=1, layer=0, complicate=0, prints=None):
        self.operator = operator
        self.layer = layer
        self.complicate = complicate
        self.isprint = isprint
        self.sub = [Sub(l, isprint, s, layer, complicate) for l, s in zip(sub, status)]
        self.values = self.get_values()
        self.prints = self.get_prints(prints)

    def get_values(self):
        for s in self.sub:
            if not s.values:
                return []

        def _get_concate(subs):
            if len(subs) == 1:
                return [[v] for v in subs[0].values]
            return [[v]+l for v in subs[0].values for l in _get_concate(subs[1:])]

        return [self._get_val(v) for v in _get_concate(self.sub)]

    def get_prints(self, _prints):
        if not self.isprint:
            return None
        for s in self.sub:
            if not s.values:
                return None
        def _get_concate(subs):
            if len(subs) == 1:
                return subs[0].prints
            return [p+l for p in subs[0].prints for l in _get_concate(subs[1:])]
        return [_prints+p for p in _get_concate(self.sub)]

    def calculate(self, x1, x2, operator):
        if operator == '+':
            return x1 + x2
        if operator == '*':
            return x1 * x2
        if operator == '/':
            return x1 / x2
        if operator == '-':
            return x1 - x2
        if operator == '**':
            return x1 ** x2
        if operator == '&1':
            return x1
        if operator == '&2':
            return x2

    def _get_val(self, data):
        subject = self.operator.split('.')
        x1 = data[int(subject[0])-1]
        x2 = data[int(subject[2])-1]
        operator1 = subject[1]
        if len(subject) == 3:
            return self.calculate(x1, x2, operator1)
        else:
            x3 = data[int(subject[4])-1]
            operator2 = subject[3]
            if operator2 in ['/', '*', '**']:
                return self.calculate(x1, self.calculate(x2, x3, operator2), operator1)
            else:
                return self.calculate(self.calculate(x1, x2, operator1), x3, operator2)


class Sub(object):

    def __init__(self, l, isprint=1, status=1, layer=1, complicate=1):
        self.l = l
        self.l_length = len(self.l)
        self.isprint = isprint
        self.status = status
        self.layer = layer
        self.complicate = complicate
        self.subs = list()
        self.values = list()
        self.prints = list()
        self.get_sub_and_val()

    def get_sub_and_val(self):
        if not self.status:
            self.values = [self.l[-1]]
            self.prints.append('')
            return 0
        else:
            for m in list(filter(lambda x: x.startswith('_is'), dir(self))):
                tep = getattr(self, m)()
                if tep is not None:
                    self.values.append(tep)
                    return 0
        res = 0
        for i in self.l:
            if not isinstance(i, int):
                return 0
        for m in list(filter(lambda x: x.startswith('_get_sub'), dir(self))):
            res += getattr(self, m)()
        if self.subs:
            for s in self.subs:
                if s.values:
                    self.values.extend(s.values)
                    if self.isprint:
                        self.prints.extend(s.prints)
        return res

    def _is_single(self):
        if self.l_length < 2:
            return None
        if len(set(self.l)) == 1:
            if self.isprint:
                self.prints.append('')
            return self.l[0]
        else:
            return None

    def _is_ap(self):
        if self.l_length < 3 or len(set(self.l)) == 1:
            return None
        sub_l = [self.l[i+1]-self.l[i] for i in range(len(self.l)-1)]
        if len(set(sub_l)) == 1:
            pad = sub_l[0]
            if self.isprint:
                self.prints.append(' '.join(map(str, self.l))+'构成等差数列，下一个数是'+str(self.l[-1]+pad)+'\n')
            return self.l[-1]+pad
        else:
            return None

    def _is_gp(self):
        if self.l_length < 3 or 0 in self.l or len(set(self.l)) == 1:
            return None
        sub_l = [self.l[i+1]/self.l[i] for i in range(self.l_length-1)]
        if len(set(sub_l)) == 1:
            pad = list(sub_l)[0]
            if self.isprint:
                self.prints.append(' '.join(map(str, self.l)) + '构成等比数列，下一个数是%d \n' % (self.l[-1]*pad))
            return int(self.l[-1]*pad)
        else:
            return None

    def _is_repeat(self):
        if self.l_length < 3 or len(set(self.l)) == 1:
            return None
        for i in range(2, int((self.l_length+1)/2)+1):
            tail = self.l_length % i
            if tail == 0:
                l_np = np.array(self.l)
            else:
                l_np = np.array(self.l)[:-tail]
            tep = ["".join(map(str, r)) for r in l_np.reshape(-1, i)]
            if len(set(tep)) == 1:
                if tail and not "".join(map(str, self.l[:tail])) == "".join(map(str, self.l[-tail:])):
                    continue
                else:
                    if self.isprint:
                        self.prints.append(' '.join(map(str, self.l))+'重复循环，下一个数是'+str(self.l[tail])+'\n')
                    return self.l[tail]
        return None

    def _get_sub_one_by_divide(self):
        if self.layer > 1 or self.complicate > 1:
            return 0
        if self.l_length < 3:
            return 0
        l1 = [[]]
        l2 = [[]]
        for n in self.l:
            if n > 0:
                indexs = range(int(pow(n, 1 / 3)), int(pow(n, 1 / 2)) + 1)
            else:
                return 0
            l1 = [_v + [i] for i in indexs for _v in l1 if not n % i]
            l2 = [_v + [n // i] for i in indexs for _v in l2 if not n % i]
        for _ in range(len(l1)):
            prints = ''
            if self.isprint:
                for i in range(self.l_length):
                    prints = prints + str(l1[_][i]) + 'x' + str(l2[_][i]) + '=' + str(self.l[i]) + '; '
                prints = prints + '\n'
            self.subs.append(Sequence(operator='1.*.2', sub=[l1[_], l2[_]], status=[1, 1], isprint=self.isprint,
                                      layer=self.layer + 1, complicate=self.complicate + 1, prints=prints))
        return 1

    def _get_sub_one_by_sqrt(self):
        if self.layer > 2 or self.complicate > 1:
            return 0
        if self.l_length < 3 or min(self.l) < 0 or max(self.l) < 4:
            return 0
        l1 = [[]]
        l2 = [[]]
        for i, n in enumerate(self.l):
            if n == -1:
                l1 = [_v + [0] for _v in l1]
                l2 = [_v + [-1] for _v in l2]
            else:
                candi = []
                tep = int(pow(n, 1 / 2))
                dif_l = n - pow(tep, 2)
                if dif_l == 0:
                    candi.append([tep, dif_l])
                else:
                    if dif_l <= i + 2:
                        candi.append([tep, dif_l])
                    dif_u = pow(tep + 1, 2) - n
                    if dif_u <= i + 2:
                        candi.append([tep + 1, -dif_u])
                if candi:
                    l1 = [_v + [c[0]] for c in candi for _v in l1]
                    l2 = [_v + [c[1]] for c in candi for _v in l2]
                else:
                    return 0
        for _ in range(len(l1)):
            prints = ''
            if self.isprint:
                for i in range(self.l_length):
                    if l2[_][i] >= 0:
                        add_ = '+'
                    else:
                        add_ = ''
                    prints = prints + str(l1[_][i])+'^2'+add_+str(l2[_][i])+'='+str(self.l[i])+'; '
                prints = prints + '\n'
            self.subs.append(Sequence(operator='1.**.2.+.3', sub=[l1[_], [2]*self.l_length, l2[_]],
                                      status=[1, 0, 1], isprint=self.isprint,
                                      layer=self.layer+1, complicate=self.complicate+1, prints=prints))
        return 1

    def _get_sub_one_by_cube_root(self):
        if self.layer > 2 or self.complicate > 1:
            return 0
        if self.l_length < 3 or max(self.l) < 4:
            return 0
        l1 = [[]]
        l2 = [[]]
        for i, n in enumerate(self.l):
            if n == 0:
                l1 = [_v + [c] for c in [-1, 0] for _v in l1]
                l2 = [_v + [c] for c in [1, 0] for _v in l2]
                continue
            if n == 1:
                l1 = [_v + [c] for c in [1, 0, -1] for _v in l1]
                l2 = [_v + [c] for c in [-1, 0, 1] for _v in l2]
                continue
            if n == -1:
                l1 = [_v + [c] for c in [0, -1] for _v in l1]
                l2 = [_v + [c] for c in [-1, 0] for _v in l2]
            if n < 0:
                sigm = -1
            else:
                sigm = 1
            candi = []
            tep = int(pow(n*sigm, 1 / 3)+0.499*(1-sigm))*sigm
            dif_l = n - pow(tep, 3)
            if dif_l == 0:
                candi.append([tep, dif_l])
            else:
                if dif_l <= i*2 + 3:
                    candi.append([tep, dif_l])
                dif_u = pow(tep + 1, 3) - n
                if dif_u <= i*2 + 3:
                    candi.append([tep + 1, -dif_u])
            if candi:
                l1 = [_v + [c[0]] for c in candi for _v in l1]
                l2 = [_v + [c[1]] for c in candi for _v in l2]
            else:
                return 0
        for _ in range(len(l1)):
            prints = ''
            if self.isprint:
                for i in range(self.l_length):
                    if l2[_][i] >= 0:
                        add_ = '+'
                    else:
                        add_ = ''
                    prints = prints + str(l1[_][i])+'^3'+add_+str(l2[_][i])+'='+str(self.l[i])+'; '
                prints = prints + '\n'
            self.subs.append(Sequence(operator='1.**.2.+.3', sub=[l1[_], [3]*self.l_length, l2[_]],
                                      status=[1, 0, 1], isprint=self.isprint,
                                      layer=self.layer+1, complicate=self.complicate+1, prints=prints))
        return 1

    def __get_sub_one_by_gcd(self):
        if self.layer > 1:
            return 0

        def gcd(a,b):
            if a == 1 or b == 1:
                return 1
            while not a == b:
                if a > b:
                    a = a - b
                else:
                    b = b - a
            return a

        res = reduce(gcd, self.l)
        if res == 1:
            return 0
        else:
            c_list = [n//res for n in self.l]
            if self.isprint:
                prints = '提取最大公约数，得到'+' '.join(map(str, c_list))+'\n'
            else:
                prints = ''
            self.subs.append(Sequence(operator='1.*.2', sub=[c_list, [res]*self.l_length],
                                      status=[1, 0], isprint=self.isprint,
                                      layer=self.layer + 1, complicate=self.complicate, prints=prints))
            return 1


    def _get_sub_two_by_difference(self):
        if self.l_length < 4 or self.layer > 4:
            return 0
        c_list = [self.l[i+1]-self.l[i] for i in range(self.l_length-1)]
        if self.isprint:
            prints = ' '.join(map(str, self.l))+'，后一项减前一项，得到'+' '.join(map(str, c_list))+'\n'
        else:
            prints = ''
        self.subs.append(Sequence(operator='1.+.2', sub=[c_list, self.l], status=[1, 0], isprint=self.isprint,
                                  layer=self.layer+1, complicate=self.complicate, prints=prints))
        return 1

    def _get_sub_two_by_difference_1(self):
        if self.l_length < 4 or self.complicate > 1:
            return 0
        for t in range(2, max(3, 5-self.layer)):
            c_list = [self.l[i+1]-self.l[i]*t for i in range(self.l_length-1)]
            if self.isprint:
                prints = ' '.join(map(str, self.l))+'，后一项减前一项的'+str(t)+'倍，得到'+' '.join(map(str, c_list))+'\n'
            else:
                prints = ''
            self.subs.append(Sequence(operator='1.+.2.*.3', sub=[c_list, [2]*self.l_length, self.l], status=[1, 0, 0],
                                      isprint=self.isprint,
                                      layer=self.layer+1, complicate=self.complicate+1, prints=prints))
        return 1

    def _get_sub_two_by_difference_2(self):
        if self.l_length < 4 or self.layer > 1:
            return 0
        c_list = [self.l[i+1]-self.l[i]*(i+1) for i in range(self.l_length-1)]
        prints = ''
        if self.isprint:
            for i in range(self.l_length-1):
                prints = prints + str(self.l[i+1])+'-'+str(self.l[i])+'*'+str(i+1)+'='+str(c_list[i]) + '; '
            prints = prints + '\n'
        self.subs.append(Sequence(operator='1.+.2.*.3', sub=[c_list, self.l, list(range(1, self.l_length+1))],
                                  status=[1, 0, 0], isprint=self.isprint,
                                  layer=self.layer+1, complicate=self.complicate+1, prints=prints))
        return 1

    def _get_sub_two_by_difference_pow(self):
        if self.l_length < 4 or self.layer > 1 or self.complicate > 1:
            return 0
        if max(self.l) > 10:
            indexs = [2]
        else:
            indexs = [2,3]
        for j in indexs:
            c_list = [self.l[i+1]-self.l[i]**j for i in range(self.l_length-1)]
            if self.isprint:
                prints = ' '.join(map(str, self.l))+'，后一项减前一项的'+str(j)+'次方，得到'+' '.join(map(str, c_list)) + '\n'
            else:
                prints = ''
            self.subs.append(Sequence(operator='1.+.2.**.3', sub=[c_list, self.l, [j]*self.l_length], status=[1, 0, 0],
                                      isprint=self.isprint,
                                      layer=self.layer+1, complicate=self.complicate+1, prints=prints))
        return 1

    def _get_sub_two_by_add(self):
        if self.layer > 1:
            return 0
        if self.l_length < 4:
            return 0
        c_list = [self.l[i]+self.l[i+1] for i in range(self.l_length-1)]
        if self.isprint:
            prints = '相邻两项相加，得到'+' '.join(map(str, c_list)) + '\n'
        else:
            prints = ''
        self.subs.append(Sequence(operator='1.-.2', sub=[c_list, self.l], status=[1, 0], isprint=self.isprint,
                                  layer=self.layer+1, complicate=self.complicate+1, prints=prints))
        return 1

    def _get_sub_two_by_division(self):
        if self.layer > 1:
            return 0
        if self.l_length < 4 or 0 in self.l:
            return 0
        c_list = []
        for i in range(self.l_length-1):
            tep = self.l[i+1] / self.l[i]
            if int(tep*10)-tep*10 == 0:
                c_list.append(tep)
            else:
                return 0
        if self.isprint:
            prints = '相邻两项相除，得到'+' '.join(map(str, c_list)) + '\n'
        else:
            prints = ''
        self.subs.append(Sequence(operator='1.*.2', sub=[c_list, self.l], status=[1, 0], isprint=self.isprint,
                                  layer=self.layer+1, complicate=self.complicate+1, prints=prints))
        return 1




    def _get_sub_three_by_add(self):
        if self.layer > 1:
            return 0
        if self.l_length < 5:
            return 0
        c_list = [self.l[i]+self.l[i+1]+self.l[i+2] for i in range(len(self.l)-2)]
        if self.isprint:
            prints = '相邻三项相加，得到'+' '.join(map(str, c_list)) + '\n'
        else:
            prints = ''
        self.subs.append(Sequence(operator='1.-.2.-.3', sub=[c_list, self.l[1:], self.l[2:]], status=[1, 0, 0],
                                  isprint=self.isprint,
                                  layer=self.layer+1, complicate=self.complicate+1, prints=prints))
        return 1

    def _get_sub_three_by_diff_pow(self):
        if self.layer > 1 or self.complicate > 1:
            return 0
        if self.l_length < 4:
            return 0
        for i in range(1,3):
            for j in range(1,3):
                if i == j == 1:
                    continue
                l1 = [self.l[k]**i+self.l[k+1]**j for k in range(len(self.l)-1)]
                l2 = [self.l[k+2]-l1[k] for k in range(len(l1)-1)]
                if self.isprint:
                    prints = '第三项减去第一项的%s次方再减去第二项的%s次方，得到' % (i, j) + ' '.join(map(str, l2)) + '\n'
                else:
                    prints = ''
                self.subs.append(Sequence(operator='1.+.2', sub=[l1, l2], status=[0, 1], isprint=self.isprint,
                                          layer=self.layer+1, complicate=self.complicate+1, prints=prints))
        return 1

    def _get_sub_three_by_diff_pow2(self):
        if self.layer > 1 or self.complicate > 1:
            return 0
        if self.l_length < 4:
            return 0
        for i in range(1,3):
            for j in range(1,3):
                if i == j == 1:
                    continue
                l1 = [self.l[k]**i-self.l[k+1]**j for k in range(len(self.l)-1)]
                l2 = [self.l[k+2]-l1[k] for k in range(len(l1)-1)]
                if self.isprint:
                    prints = '第三项减去第一项的%s次方再加上第二项的%s次方，得到' % (i, j) + ' '.join(map(str, l2)) + '\n'
                else:
                    prints = ''
                self.subs.append(Sequence(operator='1.+.2', sub=[l1, l2], status=[0, 1], isprint=self.isprint,
                                          layer=self.layer+1, complicate=self.complicate+1, prints=prints))
        return 1

    def _get_sub_three_by_diff(self):
        if self.layer > 1 or self.complicate > 1:
            return 0
        if self.l_length < 4:
            return 0
        for i in range(1,4):
            for j in range(1,4):
                l1 = [self.l[k]*i+self.l[k+1]*j for k in range(len(self.l)-1)]
                l2 = [self.l[k+2]-l1[k] for k in range(len(l1)-1)]
                if self.isprint:
                    prints = '第三项减去第一项的%s倍再减去第二项的%s倍，得到' % (i, j) + ' '.join(map(str, l2)) + '\n'
                else:
                    prints = ''
                self.subs.append(Sequence(operator='1.+.2', sub=[l1, l2], status=[0, 1], isprint=self.isprint,
                                          layer=self.layer+1, complicate=self.complicate+1, prints=prints))
        return 1

    def _get_sub_four_by_diff(self):
        if self.layer > 1 or self.complicate > 1:
            return 0
        if self.l_length < 5:
            return 0
        l1 = [self.l[k+2]+self.l[k+1]+self.l[k] for k in range(self.l_length-2)]
        l2 = [self.l[k+3]-l1[k] for k in range(len(l1)-1)]
        if not len(set(l2)) == 1:
            return 0
        if self.isprint:
            prints = '前三项的和等于第四项' + '\n'
        else:
            prints = ''
        self.subs.append(Sequence(operator='1.+.2', sub=[l1, l2], status=[0, 1], isprint=self.isprint,
                                  layer=self.layer + 1, complicate=self.complicate+1, prints=prints))
        return 1


    def _get_sub_by_odd_and_even(self):
        if self.layer > 2:
            return 0
        if self.l_length < 4:
            return 0
        l1 = self.l[::2]
        l2 = self.l[1::2]
        if len(l1) == len(l2):
            operator = '1.&1.2'
        else:
            operator = '1.&2.2'
        if self.isprint:
            prints = ' '.join(map(str, self.l))+'奇数列和偶数列拆开，得到'+' '.join(map(str, l1))+'和'+' '.join(map(str, l2)) + '\n'
        else:
            prints = ''
        self.subs.append(Sequence(operator=operator, sub=[l1, l2], status=[1, 1], isprint=self.isprint,
                                  layer=self.layer+1, complicate=self.complicate+1, prints=prints))
        return 1


def find_pattern(_s, isprint=1):
    sub_ = Sub(_s, isprint)
    if not sub_.values:
        print('找不到！')
    if isprint:
        for v, s in zip(sub_.values, sub_.prints):
            print('答案：', v)
            print('过程：')
            print(s)
            print('-'*50)
    else:
        print(set(sub_.values))

if __name__ == "__main__":
    find_pattern([1,0,1,0,1], 0)
