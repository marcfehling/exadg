/*
 * template_precursor.h
 *
 */

#ifndef APPLICATIONS_INCOMPRESSIBLE_NAVIER_STOKES_TEST_CASES_TEMPLATE_PRECURSOR_H_
#define APPLICATIONS_INCOMPRESSIBLE_NAVIER_STOKES_TEST_CASES_TEMPLATE_PRECURSOR_H_

namespace IncNS
{
namespace TemplatePrecursor
{
template<int dim, typename Number>
class Application : public ApplicationBasePrecursor<dim, Number>
{
public:
  Application() : ApplicationBasePrecursor<dim, Number>("")
  {
  }

  Application(std::string input_file) : ApplicationBasePrecursor<dim, Number>(input_file)
  {
    // parse application-specific parameters
    ParameterHandler prm;
    this->add_parameters(prm);
    parse_input(input_file, prm, true, true);
  }

  void
  add_parameters(ParameterHandler & prm)
  {
    // clang-format off
    prm.enter_subsection("Application");
      // TODO
    prm.leave_subsection();
    // clang-format on
  }

  void
  set_input_parameters(InputParameters & param)
  {
    (void)param;
  }

  void
  set_input_parameters_precursor(InputParameters & param)
  {
    (void)param;
  }

  void
  create_grid(std::shared_ptr<parallel::TriangulationBase<dim>> triangulation,
              unsigned int const                                n_refine_space,
              std::vector<GridTools::PeriodicFacePair<typename Triangulation<dim>::cell_iterator>> &
                periodic_faces)
  {
    (void)triangulation;
    (void)n_refine_space;
    (void)periodic_faces;
  }

  void
  create_grid_precursor(
    std::shared_ptr<parallel::TriangulationBase<dim>> triangulation,
    unsigned int const                                n_refine_space,
    std::vector<GridTools::PeriodicFacePair<typename Triangulation<dim>::cell_iterator>> &
      periodic_faces)
  {
    (void)triangulation;
    (void)n_refine_space;
    (void)periodic_faces;
  }

  void
  set_boundary_conditions(std::shared_ptr<BoundaryDescriptorU<dim>> boundary_descriptor_velocity,
                          std::shared_ptr<BoundaryDescriptorP<dim>> boundary_descriptor_pressure)
  {
    typedef typename std::pair<types::boundary_id, std::shared_ptr<Function<dim>>> pair;

    (void)boundary_descriptor_velocity;
    (void)boundary_descriptor_pressure;
  }

  void
  set_boundary_conditions_precursor(
    std::shared_ptr<BoundaryDescriptorU<dim>> boundary_descriptor_velocity,
    std::shared_ptr<BoundaryDescriptorP<dim>> boundary_descriptor_pressure)
  {
    typedef typename std::pair<types::boundary_id, std::shared_ptr<Function<dim>>> pair;

    (void)boundary_descriptor_velocity;
    (void)boundary_descriptor_pressure;
  }

  void
  set_field_functions(std::shared_ptr<FieldFunctions<dim>> field_functions)
  {
    field_functions->initial_solution_velocity.reset(new Functions::ZeroFunction<dim>(dim));
    field_functions->initial_solution_pressure.reset(new Functions::ZeroFunction<dim>(1));
    field_functions->analytical_solution_pressure.reset(new Functions::ZeroFunction<dim>(1));
    field_functions->right_hand_side.reset(new Functions::ZeroFunction<dim>(dim));
  }

  void
  set_field_functions_precursor(std::shared_ptr<FieldFunctions<dim>> field_functions)
  {
    field_functions->initial_solution_velocity.reset(new Functions::ZeroFunction<dim>(dim));
    field_functions->initial_solution_pressure.reset(new Functions::ZeroFunction<dim>(1));
    field_functions->analytical_solution_pressure.reset(new Functions::ZeroFunction<dim>(1));
    field_functions->right_hand_side.reset(new Functions::ZeroFunction<dim>(dim));
  }

  std::shared_ptr<PostProcessorBase<dim, Number>>
  construct_postprocessor(unsigned int const degree, MPI_Comm const & mpi_comm)
  {
    (void)degree;

    std::shared_ptr<PostProcessorBase<dim, Number>> pp;

    PostProcessorData<dim> pp_data;
    pp.reset(new PostProcessor<dim, Number>(pp_data, mpi_comm));

    return pp;
  }

  std::shared_ptr<PostProcessorBase<dim, Number>>
  construct_postprocessor_precursor(unsigned int const degree, MPI_Comm const & mpi_comm)
  {
    (void)degree;

    std::shared_ptr<PostProcessorBase<dim, Number>> pp;

    PostProcessorData<dim> pp_data;
    pp.reset(new PostProcessor<dim, Number>(pp_data, mpi_comm));

    return pp;
  }
};

} // namespace TemplatePrecursor
} // namespace IncNS

#endif /* APPLICATIONS_INCOMPRESSIBLE_NAVIER_STOKES_TEST_CASES_TEMPLATE_PRECURSOR_H_ */
