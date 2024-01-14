오후 1:05 2023-11-29

<수정한 내역들>
1. Feasibility check
- Feasibility problem에서 도출하는 W가 rank-1을 만족하지 않는 것을 확인
	(FP1에서 W는 rank-1 만족해야함)
- 논문에서 주장하는 대로 수정 
	(W를 0으로 initialize하고 constraint 만족하는 R 찾는 식으로)

2. rank-1 optimality cond.
- CVX로부터 도출한 솔루션을 rank-1으로 바꿔주는 과정 추가
- 논문이 주장하는 과정에서 feasibility 유지되고 optimal function 값 유지되는 것 확인

3. 탈출조건 변경
- n번째 episode에서 \hat{r_k}의 값 변화량이 1e-6보다 작을 때 탈출하도록 수정


<변화>
- 일단 NaN 뜨는 부분들은 수정
- 한 쪽으로 솔루션이 계속 몰리는데, 이 부분은 매 iteration마다 뜯어보니까 항상 한 쪽으로 몰림

- Comm. only보다 radar SNR const.가 -13dB일 때 rate가 더 높게 나옴 <- 이 부분이 아직 이상하다