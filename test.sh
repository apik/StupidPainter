NC=3

echo -e \
     "[T(1)T(1)]\n" \
     "f(1,2,3)f(1,2,3)\n" \
     "[T(1)T(1)T(2)T(2)]\n" \
     "[T(1)T(2)T(1)T(2)]\n" \
     "[T(1)T(2)][T(1)T(2)]\n" \
     "[T(1)T(2)T(3)][T(1)T(2)T(3)]\n" \
    | ./StupidPainter $NC
