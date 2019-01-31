#include "constraint.h"

TestConstraint::TestConstraint(int id_0, int id_1, int id_2, const VectorX &x, ScalarType rangemin, ScalarType rangemax, ScalarType weight):
	Constraint(CONSTRAINT_TYPE_TRI) {//initial set for constrain
									 //ids 0,1,2
	range_min = rangemin;
	range_max = rangemax;
	EigenMatrix32 edge, P;

	ids[0] = id_0;
	ids[1] = id_1;
	ids[2] = id_2;
	EigenVector3 pos_0, pos_1, pos_2;
	pos_0 = { x[3 * ids[0]] ,x[3 * ids[0] + 1] ,x[3 * ids[0] + 2] };
	pos_1 = { x[3 * ids[1]] ,x[3 * ids[1] + 1] ,x[3 * ids[1] + 2] };
	pos_2 = { x[3 * ids[2]] ,x[3 * ids[2] + 1] ,x[3 * ids[2] + 2] };
	edge.col(0) = pos_1 - pos_0;
	edge.col(1) = pos_2 - pos_0;

	P.col(0) = edge.col(0).normalized();
	P.col(1) = (edge.col(1) - edge.col(1).dot(P.col(0))*P.col(0)).normalized();

	rest = (P.transpose()*edge).inverse();

	Eigen::Matrix<ScalarType, 2, 3> IND;
	IND.block<2, 1>(0, 0) = EigenVector2(-1, -1);
	IND.block<2, 2>(0, 1) = EigenMatrix2::Identity();
	m_G = rest.transpose() * IND;
	ScalarType Aear = ((P.transpose()*edge)).determinant() / 2.0;
	weight *= std::sqrt(std::abs(Aear));
	//set the normal vector 

}
TestConstraint::TestConstraint(int id_0, int id_1, int id_2, const std::vector<float> &x, ScalarType rangemin, ScalarType rangemax, ScalarType weight) :
	Constraint(CONSTRAINT_TYPE_TRI) {//initial set for constrain
									 //ids 0,1,2
	range_min = rangemin;
	range_max = rangemax;
	EigenMatrix32 edge, P;

	ids[0] = id_0;
	ids[1] = id_1;
	ids[2] = id_2;
	EigenVector3 pos_0, pos_1, pos_2;
	pos_0 = { x[3 * ids[0]] ,x[3 * ids[0] + 1] ,x[3 * ids[0] + 2] };
	pos_1 = { x[3 * ids[1]] ,x[3 * ids[1] + 1] ,x[3 * ids[1] + 2] };
	pos_2 = { x[3 * ids[2]] ,x[3 * ids[2] + 1] ,x[3 * ids[2] + 2] };
	edge.col(0) = pos_1 - pos_0;
	edge.col(1) = pos_2 - pos_0;

	P.col(0) = edge.col(0).normalized();
	P.col(1) = (edge.col(1) - edge.col(1).dot(P.col(0))*P.col(0)).normalized();

	rest = (P.transpose()*edge).inverse();

	Eigen::Matrix<ScalarType, 2, 3> IND;
	IND.block<2, 1>(0, 0) = EigenVector2(-1, -1);
	IND.block<2, 2>(0, 1) = EigenMatrix2::Identity();
	m_G = rest.transpose() * IND;
	ScalarType Aear = ((P.transpose()*edge)).determinant() / 2.0;
	weight *= std::sqrt(std::abs(Aear));
}

TestConstraint::TestConstraint(int id_0, int id_1, int id_2, const VectorX &_temp_uv_added,ScalarType rangemin, ScalarType rangemax, ScalarType weight,int tem) :
	Constraint(CONSTRAINT_TYPE_TRI) {//initial set for constrain
									 //ids 0,1,234e5
	range_min = rangemin;
	range_max = rangemax;
	EigenMatrix32 edge, P;

	ids[0] = id_0;
	ids[1] = id_1;
	ids[2] = id_2;
	EigenVector3 pos_0, pos_1, pos_2;
	pos_0 = { _temp_uv_added[3 * ids[0]] ,_temp_uv_added[3 * ids[0] + 1] ,_temp_uv_added[3 * ids[0] + 2] };
	pos_1 = { _temp_uv_added[3 * ids[1]] ,_temp_uv_added[3 * ids[1] + 1] ,_temp_uv_added[3 * ids[1] + 2] };
	pos_2 = { _temp_uv_added[3 * ids[2]] ,_temp_uv_added[3 * ids[2] + 1] ,_temp_uv_added[3 * ids[2] + 2] };

	edge.col(0) = pos_1 - pos_0;
	edge.col(1) = pos_2 - pos_0;

	P.col(0) = edge.col(0).normalized();
	P.col(1) = (edge.col(1) - edge.col(1).dot(P.col(0))*P.col(0)).normalized();

	rest = (P.transpose()*edge).inverse();
	//Eigen::Matrix<ScalarType, 2, 3> IND;
	//IND.block<2, 1>(0, 0) = EigenVector2(-1, -1);
	//IND.block<2, 2>(0, 1) = EigenMatrix2::Identity();
	//m_G = rest.transpose() * IND;
	ScalarType Aear = ((P.transpose()*edge)).determinant() / 2.0;
	weight *= std::sqrt(std::abs(Aear));
}

void TestConstraint::projection_uvs(const VectorX &x, EigenMatrix2X & projection, int id)//have to add the ides for 3indez

{
	EigenMatrix32 edge, P;

	EigenVector3 pos_0, pos_1, pos_2;
	pos_0 = { x[3 * ids[0]] ,x[3 * ids[0] + 1] ,x[3 * ids[0] + 2] };
	pos_1 = { x[3 * ids[1]] ,x[3 * ids[1] + 1] ,x[3 * ids[1] + 2] };
	pos_2 = { x[3 * ids[2]] ,x[3 * ids[2] + 1] ,x[3 * ids[2] + 2] };
	edge.col(0) = pos_1 - pos_0;
	edge.col(1) = pos_2 - pos_0;


	P.col(0) = edge.col(0).normalized();
	P.col(1) = (edge.col(1) - edge.col(1).dot(P.col(0))*P.col(0)).normalized();
	EigenMatrix22 F = P.transpose()* edge* rest;

	//Dm_inv = edge;
	Eigen::JacobiSVD<EigenMatrix22>  svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);//
	EigenVector2 S = svd.singularValues();

	S(0) = clamp(S(0), range_min, range_max);
	S(1) = clamp(S(1), range_min, range_max);

	F = svd.matrixU()*S.asDiagonal()*svd.matrixV().transpose();
	projection.block<2, 2>(0, id * 2) = (F);


}


void TestConstraint::projection_set(const VectorX &x, EigenMatrix2X & projection, int id)//have to add the ides for 3indez

{
	EigenMatrix32 edge, P;

	EigenVector3 pos_0, pos_1, pos_2;
	pos_0 = { x[3 * ids[0]] ,x[3 * ids[0] + 1] ,x[3 * ids[0] + 2] };
	pos_1 = { x[3 * ids[1]] ,x[3 * ids[1] + 1] ,x[3 * ids[1] + 2] };
	pos_2 = { x[3 * ids[2]] ,x[3 * ids[2] + 1] ,x[3 * ids[2] + 2] };

	edge.col(0) = pos_1 - pos_0;
	edge.col(1) = pos_2 - pos_0;

	P.col(0) = edge.col(0).normalized();
	P.col(1) = (edge.col(1) - edge.col(1).dot(P.col(0))*P.col(0)).normalized();
	EigenMatrix22 F = P.transpose()* edge* rest;

	Eigen::JacobiSVD<EigenMatrix22>  svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);//
	EigenVector2 S = svd.singularValues();

	S(0) = clamp(S(0), range_min, range_max);
	S(1) = clamp(S(1), range_min, range_max);
	
	F = svd.matrixU()*S.asDiagonal()*svd.matrixV().transpose();//F(22) = svd (22) * S(norm*svd(22)

	projection.block<2, 2>(0, id * 2) = (F);//P (3*2)*(2*") = (3*2)

}

void TestConstraint::projection_set_after(const VectorX &x, EigenMatrix2X & projection, EigenMatrix3X & project, int id)//
{
	EigenMatrix32 edge, P;
	EigenVector3 pos_0, pos_1, pos_2;
	pos_0 = { x[3 * ids[0]] ,x[3 * ids[0] + 1] ,x[3 * ids[0] + 2] };
	pos_1 = { x[3 * ids[1]] ,x[3 * ids[1] + 1] ,x[3 * ids[1] + 2] };
	pos_2 = { x[3 * ids[2]] ,x[3 * ids[2] + 1] ,x[3 * ids[2] + 2] };

	edge.col(0) = pos_1 - pos_0;
	edge.col(1) = pos_2 - pos_0;

	P.col(0) = edge.col(0).normalized();
	P.col(1) = (edge.col(1) - edge.col(1).dot(P.col(0))*P.col(0)).normalized();
	EigenMatrix22 F = projection.block<2, 2>(0, id * 2);
	project.block<3, 2>(0, id * 2) = (P*F);//P (3*2)*(2*") = (3*2)
}

void TestConstraint::projection(const VectorX &x, EigenMatrix3X & project, int id)//have to add the ides for 3indez
{
	EigenMatrix32 edge, P;

	EigenVector3 pos_0, pos_1, pos_2;
	pos_0 = { x[3 * ids[0]] ,x[3 * ids[0] + 1] ,x[3 * ids[0] + 2] };
	pos_1 = { x[3 * ids[1]] ,x[3 * ids[1] + 1] ,x[3 * ids[1] + 2] };
	pos_2 = { x[3 * ids[2]] ,x[3 * ids[2] + 1] ,x[3 * ids[2] + 2] };

	edge.col(0) = pos_1 - pos_0;
	edge.col(1) = pos_2 - pos_0;

	P.col(0) = edge.col(0).normalized();
	P.col(1) = (edge.col(1) - edge.col(1).dot(P.col(0))*P.col(0)).normalized();
	EigenMatrix22 F = P.transpose()* edge* rest;

	Eigen::JacobiSVD<EigenMatrix22>  svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);//
	EigenVector2 S = svd.singularValues();

	S(0) = clamp(S(0), range_min, range_max);
	S(1) = clamp(S(1), range_min, range_max);

	F = svd.matrixU()*S.asDiagonal()*svd.matrixV().transpose();//F(22) = svd (22) * S(norm*svd(22)

	project.block<3, 2>(0, id * 2) = (P*F);//P (3*2)*(2*") = (3*2)

}


void TestConstraint::EvaluateWeightedLaplacian1D(std::vector<SparseMatrixTriplet>& laplacian_1d_triplets)
{
	ScalarType ks = 1.0;// stiffness;// m_laplacian_coeff*m_w;
	EigenMatrix3 L = ks * m_G.transpose()*m_G;

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			laplacian_1d_triplets.push_back(SparseMatrixTriplet(ids[i], ids[j], L(i, j)));
		}
	}
	//set the laplacian (after changed)
}


void TestConstraint::EvaluateJMatrix(unsigned int index, std::vector<SparseMatrixTriplet>& J_triplets)
{
	// p is th base	
	J_triplets.push_back(SparseMatrixTriplet(2 * index + 0, ids[0], m_G(0, 0)));
	J_triplets.push_back(SparseMatrixTriplet(2 * index + 0, ids[1], m_G(0, 1)));
	J_triplets.push_back(SparseMatrixTriplet(2 * index + 1, ids[1], m_G(1, 1)));
	J_triplets.push_back(SparseMatrixTriplet(2 * index + 1, ids[0], m_G(1, 0)));
	J_triplets.push_back(SparseMatrixTriplet(2 * index + 0, ids[2], m_G(0, 2)));
	J_triplets.push_back(SparseMatrixTriplet(2 * index + 1, ids[2], m_G(1, 2)));
}

double TestConstraint::clamp(double n, double lower, double upper) {
	//this is not for
	return std::max(lower, std::min(n, upper));
	//here to the normal

}


void TestConstraint::projection_interpolate_normal(const VectorX &x, EigenMatrix2X & projection, EigenMatrix3X & project, int id, std::vector<int> triangle_list)//
{
	EigenMatrix32 edge, P;
	EigenVector3 pos_0, pos_1, pos_2;
	pos_0 = { x[3 * ids[0]] ,x[3 * ids[0] + 1] ,x[3 * ids[0] + 2] };
	pos_1 = { x[3 * ids[1]] ,x[3 * ids[1] + 1] ,x[3 * ids[1] + 2] };
	pos_2 = { x[3 * ids[2]] ,x[3 * ids[2] + 1] ,x[3 * ids[2] + 2] };

	edge.col(0) = pos_1 - pos_0;
	edge.col(1) = pos_2 - pos_0;

	P.col(0) = edge.col(0).normalized();
	P.col(1) = (edge.col(1) - edge.col(1).dot(P.col(0))*P.col(0)).normalized();

	//set interpolation
	int element_num = triangle_list.size();
	EigenMatrix22 temp;
	temp.setZero();
	for (int i = 0; i < element_num; i++)
	{
		int triangle_index = triangle_list[i];
		temp = temp + projection.block<2, 2>(0, triangle_index * 2);
	}
	temp = temp / double(element_num);

	EigenMatrix22 F = temp;// projection.block<2, 2>(0, id * 2);
	project.block<3, 2>(0, id * 2) = (P*F);//P (3*2)*(2*") = (3*2)
}

